"""Staged ML-potential MD (ASE) for the OxNi models.

  python run_md.py ../structures/S1_ni211_meo_pore_gas.extxyz \
      --calc mace:medium-omat-0 --device cuda --protocol ../protocols/mdr_873K.json \
      --out ../runs/S1_873K [--plumed ../plumed/ch4_S1.dat]

Protocol JSON: {"log_every": 50, "traj_every": 200, "stages": [ {...}, ... ]}
  MD stage : name, T, [T_end (linear ramp)], steps, dt_fs, thermostat (langevin|bussi|nve),
             friction_per_fs (langevin) or taut_fs (bussi), freeze_tags (list of tags)
  opt stage: {"type": "opt", "fmax": 0.2, "steps": 300, "freeze_tags": [...]}
Tags: 0 Ni, 1 Ni frozen, 2 MEO, 3 support, 4 support frozen, 5 gas, 6 pre-adsorbed.

Outputs: <out>.extxyz (trajectory incl. momenta), <out>_log.csv (T, E, species, coke and
Ni metrics), <out>_final.extxyz (restart), <out>_timing.json.
An upper harmonic wall at atoms.info['wall_z'] (or --wall-z) keeps the gas in the box.
"""
import argparse
import csv
import math
import json
import os
import sys
import time

import numpy as np
from ase import units
from ase.constraints import FixAtoms
from ase.io import read, write
from ase.md.bussi import Bussi
from ase.md.langevin import Langevin
try:
    from ase.md.velocitydistribution import thermalize_momenta
except ImportError:  # ASE < 3.26
    from ase.md.velocitydistribution import MaxwellBoltzmannDistribution as thermalize_momenta
from ase.md.verlet import VelocityVerlet
from ase.optimize import FIRE

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import species as sp  # noqa: E402
from calculators import get_calculator  # noqa: E402


class ZWall:
    """Harmonic upper wall E = k/2 (z - z0)^2 for every atom with z > z0."""

    def __init__(self, z0, k=10.0):
        self.z0, self.k = float(z0), float(k)

    def adjust_positions(self, atoms, new):
        pass

    def adjust_momenta(self, atoms, momenta):
        pass

    def adjust_forces(self, atoms, forces):
        dz = atoms.positions[:, 2] - self.z0
        m = dz > 0
        forces[m, 2] -= self.k * dz[m]

    def adjust_potential_energy(self, atoms):
        dz = atoms.positions[:, 2] - self.z0
        return 0.5 * self.k * float(np.sum(dz[dz > 0] ** 2))

    def get_removed_dof(self, atoms):
        return 0

    def index_shuffle(self, atoms, ind):
        pass

    def copy(self):
        return ZWall(self.z0, self.k)

    def todict(self):
        return {"name": "ZWall", "kwargs": {"z0": self.z0, "k": self.k}}


def keep_info(info):
    """Scalars and numeric tuples (pore centre, wall) survive into the written frames."""
    out = {}
    for k, v in info.items():
        if isinstance(v, (int, float, str)):
            out[k] = v
        elif isinstance(v, (tuple, list, np.ndarray)) and np.asarray(v).dtype.kind in "if":
            out[k] = np.asarray(v, float)
    return out


def make_dyn(atoms, st, seed):
    dt = st.get("dt_fs", 0.5) * units.fs
    th = st.get("thermostat", "langevin")
    rng = np.random.default_rng(seed)
    if th == "langevin":
        return Langevin(atoms, dt, temperature_K=st["T"], fixcm=False, rng=rng,
                        friction=st.get("friction_per_fs", 0.01) / units.fs)
    if th == "bussi":
        return Bussi(atoms, dt, temperature_K=st["T"], taut=st.get("taut_fs", 100.0) * units.fs,
                     rng=rng)
    if th == "nve":
        return VelocityVerlet(atoms, dt)
    raise ValueError(th)


def set_target_T(dyn, T):
    if isinstance(dyn, Langevin):
        dyn.set_temperature(temperature_K=T)
    elif isinstance(dyn, Bussi):
        dyn.temp = T * units.kB
        dyn.target_kinetic_energy = 0.5 * dyn.temp * dyn.ndof


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("structure")
    ap.add_argument("--calc", default="mace:medium")
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--dtype", default="float32")
    ap.add_argument("--protocol", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--plumed", default=None, help="PLUMED input (UNITS LENGTH=A ENERGY=eV TIME=fs)")
    ap.add_argument("--wall-z", type=float, default=None)
    ap.add_argument("--no-wall", action="store_true")
    ap.add_argument("--no-ni-metrics", action="store_true")
    ap.add_argument("--reset-velocities", action="store_true")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--threads", type=int, default=os.cpu_count())
    a = ap.parse_args()

    import torch
    torch.set_num_threads(a.threads)

    proto = json.load(open(a.protocol))
    stages = proto["stages"]
    log_every, traj_every = proto.get("log_every", 50), proto.get("traj_every", 200)
    out = os.path.abspath(a.out)
    os.makedirs(os.path.dirname(out), exist_ok=True)

    atoms = read(a.structure)  # extxyz with tags
    tags = atoms.get_tags()
    wall = None
    if not a.no_wall and (a.wall_z or atoms.info.get("wall_z")):
        wall = ZWall(a.wall_z or atoms.info["wall_z"])

    calc = get_calculator(a.calc, a.device, a.dtype)
    if a.plumed:
        from ase.calculators.plumed import Plumed
        dts = {st.get("dt_fs", 0.5) for st in stages if st.get("type", "md") == "md"}
        assert len(dts) == 1, "PLUMED needs one timestep for all MD stages"
        # PLUMED counts every force call as an MD step, so an optimiser would deposit bias
        assert all(st.get("type", "md") == "md" for st in stages), "no opt stages with --plumed"
        lines = [ln.rstrip() for ln in open(a.plumed) if ln.strip() and not ln.lstrip().startswith("#")]
        T0 = next(st["T"] for st in stages if st.get("type", "md") == "md")
        os.chdir(os.path.dirname(out))  # COLVAR/KERNELS land next to the trajectory
        calc = Plumed(calc=calc, input=lines, timestep=dts.pop() * units.fs, atoms=atoms,
                      kT=units.kB * T0, log=out + "_plumed.log")
    atoms.calc = calc

    first_md = next((st for st in stages if st.get("type", "md") == "md"), None)
    has_p = atoms.has("momenta") and np.abs(atoms.get_momenta()).sum() > 0
    if first_md and (a.reset_velocities or not has_p):
        thermalize_momenta(atoms, temperature_K=first_md["T"], rng=np.random.default_rng(a.seed))
        p = atoms.get_momenta()
        p[np.isin(tags, first_md.get("freeze_tags", []))] = 0.0
        atoms.set_momenta(p)

    fields = (["step", "time_ps", "stage", "T_target", "T", "epot"] + sp.KEY_SPECIES +
              ["n_C", "C_in_Cn", "largest_C_cluster", "C_dissolved", "O_lattice_in_C_species",
               "O_lattice_in_gas", "H_on_oxide", "Ni_Rg", "Ni_meanCN", "Ni_detached",
               "Ni_oxidic", "Ni_meo_to_metal", "other_species"])
    fh = open(out + "_log.csv", "w", newline="")
    wr = csv.DictWriter(fh, fieldnames=fields)
    wr.writeheader()
    traj_path = out + ".extxyz"
    if os.path.exists(traj_path):
        os.remove(traj_path)

    state = dict(step=0, time_fs=0.0, stage="", T_target=0.0)
    timing = []

    def log_row():
        T = atoms.get_temperature()
        if state["T_target"] > 0 and T > 3 * state["T_target"] + 500:
            raise RuntimeError("temperature blow-up (T = %.0f K) at step %d" % (T, state["step"]))
        spc, met = sp.analyse(atoms)
        row = dict(step=state["step"], time_ps=round(state["time_fs"] / 1000, 5),
                   stage=state["stage"], T_target=state["T_target"], T=round(T, 1),
                   epot=round(atoms.get_potential_energy(), 5))
        row.update({k: spc.pop(k, 0) for k in sp.KEY_SPECIES})
        row.update(met)
        if not a.no_ni_metrics:
            row.update({k: (round(v, 4) if isinstance(v, float) else v)
                        for k, v in sp.ni_metrics(atoms).items()})
        row["other_species"] = json.dumps(dict(spc)) if spc else ""
        wr.writerow(row)
        fh.flush()

    def write_frame():
        f = atoms.copy()
        f.set_constraint()
        f.info = keep_info(atoms.info)
        f.info.update(step=state["step"], time_ps=state["time_fs"] / 1000,
                      stage=state["stage"], epot=float(atoms.get_potential_energy()))
        write(traj_path, f, append=True)

    for k, st in enumerate(stages):
        cons = [FixAtoms(mask=np.isin(tags, st.get("freeze_tags", [])))]
        if wall is not None:
            cons.append(wall)
        atoms.set_constraint(cons)
        state["stage"] = st.get("name", "stage%d" % k)
        t0 = time.time()
        if st.get("type", "md") == "opt":
            opt = FIRE(atoms, logfile=out + "_opt_%s.log" % state["stage"])
            opt.run(fmax=st.get("fmax", 0.2), steps=st.get("steps", 300))
            state["T_target"] = 0.0
            log_row()
            write_frame()
            timing.append(dict(stage=state["stage"], type="opt", steps=opt.nsteps,
                               seconds=time.time() - t0))
            print("[%s] FIRE %d steps, fmax target %.2f" % (state["stage"], opt.nsteps,
                                                             st.get("fmax", 0.2)), flush=True)
            continue

        dyn = make_dyn(atoms, st, a.seed + k)
        dt_fs = st.get("dt_fs", 0.5)
        n, T0, T1 = st["steps"], st["T"], st.get("T_end", st["T"])
        base_step, base_time = state["step"], state["time_fs"]
        state["T_target"] = T0

        def observer():
            state["step"] = base_step + dyn.nsteps
            state["time_fs"] = base_time + dyn.nsteps * dt_fs
            if T1 != T0:
                state["T_target"] = T0 + (T1 - T0) * dyn.nsteps / n
                set_target_T(dyn, state["T_target"])
            if dyn.nsteps % log_every == 0 and (dyn.nsteps > 0 or k == 0):
                log_row()
            if dyn.nsteps % traj_every == 0 and (dyn.nsteps > 0 or k == 0):
                write_frame()

        dyn.attach(observer, interval=1 if T1 != T0 else math.gcd(log_every, traj_every))
        dyn.run(n)
        sec = time.time() - t0
        ns_day = n * dt_fs * 1e-6 / (sec / 86400)
        timing.append(dict(stage=state["stage"], steps=n, dt_fs=dt_fs, seconds=sec,
                           s_per_step=sec / n, ns_per_day=ns_day, n_atoms=len(atoms)))
        print("[%s] %d steps, %.3f s/step, %.3f ns/day, T=%.0f K" % (
            state["stage"], n, sec / n, ns_day, atoms.get_temperature()), flush=True)

    f = atoms.copy()
    f.set_constraint()
    f.info = keep_info(atoms.info)
    write(out + "_final.extxyz", f)
    json.dump(dict(calc=a.calc, device=a.device, dtype=a.dtype, stages=timing),
              open(out + "_timing.json", "w"), indent=1)
    fh.close()


if __name__ == "__main__":
    main()
