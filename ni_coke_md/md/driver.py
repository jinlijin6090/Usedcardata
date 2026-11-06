#!/usr/bin/env python
"""ReaxFF deposition / annealing driver with checkpoint-restart.

One case = one directory holding config.json.  The driver
  1. equilibrates the starting structure,
  2. inserts atoms (C and/or O) one by one at a fixed interval, aimed at the cluster,
  3. optionally deletes fragments that have left the cluster (CO, CO2, O2, stray C ...),
  4. runs a post-deposition stage (relaxation, or a long anneal when nothing is deposited).
Every `checkpoint_every_ps` it writes ckpt.data + state.json; rerunning the same command
resumes from the last checkpoint (trajectory and event log are truncated to match).

Usage:  python driver.py CASE_DIR [--threads N] [--max-hours H]
Exit codes: 0 finished, 3 stopped on --max-hours with a valid checkpoint (resubmit), other = error.
"""
import argparse
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (TYPES, ELEM, MASS, AMU_A2_FS2_TO_KCAL, KB_KCAL, fragments, formula,
                    unwrap_about_ni, write_json_atomic)

PKG_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

DEFAULTS = dict(
    T=1200.0, dt_fs=0.25, seed=1, tdamp_fs=50.0,
    equil_ps=10.0,
    deposit={"C": 0, "O": 0},
    c_interval_ps=20.0,           # time between successive C insertions; O atoms are interleaved
    post_ps=50.0,
    dump_every_ps=1.0, checkpoint_every_ps=25.0, recenter_every_ps=5.0,
    remove_desorbed=True, desorb_dist=6.0, check_desorb_every_ps=1.0,
    insert_gap=3.5, min_dist=2.2,
    ffield="ff/ffield.reax.FC",
)


def resolve(p, case_dir):
    if os.path.isabs(p):
        return p
    for base in (case_dir, PKG_ROOT):
        q = os.path.normpath(os.path.join(base, p))
        if os.path.exists(q):
            return q
    raise FileNotFoundError(
        f"{p} not found (looked relative to {case_dir} and {PKG_ROOT}). "
        "For task2_anneal this means the task1_flux case it starts from has not finished yet.")


class Case:
    def __init__(self, case_dir, threads, max_hours):
        self.dir = os.path.abspath(case_dir)
        cfg = dict(DEFAULTS)
        cfg.update(json.load(open(os.path.join(self.dir, "config.json"))))
        self.cfg = cfg
        self.threads = threads
        self.max_sec = max_hours * 3600 if max_hours else None
        self.t_wall0 = time.time()
        self.dt = cfg["dt_fs"]
        ps = lambda x: int(round(x * 1000.0 / self.dt))
        self.ps = ps

        # insertion schedule (fixed by seed, identical on every restart)
        dep = {k: int(v) for k, v in cfg["deposit"].items() if int(v) > 0}
        seq = sum(([e] * n for e, n in dep.items()), [])
        rng = np.random.default_rng(cfg["seed"])
        rng.shuffle(seq)
        self.seq = seq
        nC = dep.get("C", 0)
        span_ps = (nC if nC else len(seq)) * cfg["c_interval_ps"]
        self.ins_every = ps(span_ps / len(seq)) if seq else 0
        self.E = ps(cfg["equil_ps"])
        self.D = self.E + self.ins_every * len(seq)
        self.total = self.D + ps(cfg["post_ps"])
        self.dump_every = ps(cfg["dump_every_ps"])
        self.ckpt_every = ps(cfg["checkpoint_every_ps"])
        self.recenter_every = ps(cfg["recenter_every_ps"])
        self.desorb_every = ps(cfg["check_desorb_every_ps"])
        self.rng = np.random.default_rng(cfg["seed"] * 7919 + 17)

        self.f_state = os.path.join(self.dir, "state.json")
        self.f_traj = os.path.join(self.dir, "traj.lammpstrj")
        self.f_events = os.path.join(self.dir, "events.csv")
        self.f_ckpt = os.path.join(self.dir, "ckpt.data")
        self.f_prog = os.path.join(self.dir, "progress.log")

    # ------------------------------------------------------------------ setup
    def start_lammps(self):
        from lammps import lammps
        resumed = os.path.exists(self.f_state)
        st = json.load(open(self.f_state)) if resumed else None
        if st and st.get("phase") == "done":
            print("case already finished:", self.dir)
            sys.exit(0)
        nlog = len([f for f in os.listdir(self.dir) if f.startswith("log.lammps")])
        args = ["-log", os.path.join(self.dir, f"log.lammps.{nlog}"), "-screen", "none", "-nocite"]
        if self.threads > 1:
            args += ["-sf", "omp", "-pk", "omp", str(self.threads)]
        self.lmp = L = lammps(cmdargs=args)
        c = self.cfg
        L.command("units real")
        L.command("atom_style charge")
        L.command("boundary p p p")
        L.command("atom_modify map yes")
        if resumed:
            self.state = st
            data = self.f_ckpt
            for f, key in ((self.f_traj, "traj_bytes"), (self.f_events, "events_bytes")):
                if os.path.exists(f):
                    os.truncate(f, st[key])
        else:
            for f in (self.f_traj, self.f_events, self.f_prog):
                if os.path.exists(f):
                    os.remove(f)
            self.state = dict(step=0, j=0, phase="equil", n_restarts=0, recent={})
            data = resolve(c["init_data"], self.dir)
            with open(self.f_events, "w") as f:
                f.write("time_ps,step,event,species,nC,nO,natoms_after\n")
        L.command(f"read_data {data}")
        ff = resolve(c["ffield"], self.dir)
        L.command("pair_style reaxff NULL checkqeq yes safezone 3.0 mincap 200")
        L.command(f"pair_coeff * * {ff} C Ni O")
        L.command("fix qeq all qeq/reaxff 1 0.0 10.0 1.0e-6 reaxff")
        L.command("neighbor 2.0 bin")
        L.command("neigh_modify every 1 delay 0 check yes")
        L.command(f"timestep {self.dt}")
        L.command("compute tdyn all temp")
        L.command("compute_modify tdyn dynamic/dof yes")
        if resumed:
            L.command(f"reset_timestep {st['step']}")
            self.state["n_restarts"] = st.get("n_restarts", 0) + 1
        elif "Velocities" not in open(data).read():
            L.command(f"velocity all create {c['T']} {c['seed']*1000+3} mom yes rot yes dist gaussian")
        L.command(f"fix nvt all nvt temp {c['T']} {c['T']} {c['tdamp_fs']}")
        L.command("fix_modify nvt temp tdyn")
        L.command("fix mom all momentum 400 linear 1 1 1 angular")
        L.command("thermo_style custom step temp pe atoms")
        L.command("thermo_modify temp tdyn lost ignore norm no")
        L.command("thermo 100000000")
        L.command(f"dump d all custom {self.dump_every} {self.f_traj} id type x y z q")
        delay = st["step"] + 1 if resumed else 0
        L.command(f"dump_modify d sort id append yes delay {delay}")
        self.box = np.array([L.extract_box()[1][k] - L.extract_box()[0][k] for k in range(3)])
        self.step = self.state["step"]
        self.log(f"{'resumed at' if resumed else 'started'} step {self.step} "
                 f"(t={self.step*self.dt/1000:.2f} ps), threads={self.threads}, "
                 f"schedule: equil->{self.E}, deposit {len(self.seq)} atoms every {self.ins_every} steps "
                 f"->{self.D}, end {self.total}")

    # ------------------------------------------------------------------ helpers
    def atoms(self):
        L = self.lmp
        n = L.get_natoms()
        ids = np.array(L.numpy.extract_atom("id")[:n], dtype=int)
        types = np.array(L.numpy.extract_atom("type")[:n], dtype=int)
        pos = np.array(L.numpy.extract_atom("x")[:n], dtype=float)
        return ids, types, pos

    def log(self, msg):
        line = f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}"
        print(line, flush=True)
        with open(self.f_prog, "a") as f:
            f.write(line + "\n")

    def event(self, kind, species, nC, nO):
        with open(self.f_events, "a") as f:
            f.write(f"{self.step*self.dt/1000:.4f},{self.step},{kind},{species},{nC},{nO},"
                    f"{self.lmp.get_natoms()}\n")

    def recenter(self):
        ids, types, pos = self.atoms()
        _, c = unwrap_about_ni(pos, types, self.box)
        d = self.box / 2 - c
        self.lmp.command(f"displace_atoms all move {d[0]:.6f} {d[1]:.6f} {d[2]:.6f} units box")

    def insert(self, elem):
        ids, types, pos = self.atoms()
        pos, c = unwrap_about_ni(pos, types, self.box)
        ni = pos[types == TYPES["Ni"]]
        from scipy.spatial import cKDTree
        dni, _ = cKDTree(ni).query(pos)
        bound = dni < 3.0
        R = np.linalg.norm(pos[bound] - c, axis=1).max()
        tree = cKDTree(pos)
        r = R + self.cfg["insert_gap"]
        for attempt in range(1000):
            u = self.rng.normal(size=3)
            u /= np.linalg.norm(u)
            p = c + r * u
            if tree.query(p)[0] >= self.cfg["min_dist"]:
                break
            if attempt % 100 == 99:
                r += 0.5
        p = p - self.box * np.floor(p / self.box)           # wrap into box
        m = MASS[elem]
        sig = np.sqrt(KB_KCAL * self.cfg["T"] / (m * AMU_A2_FS2_TO_KCAL))
        speed = np.linalg.norm(self.rng.normal(scale=sig, size=3))   # Maxwell-Boltzmann speed
        v = -speed * u                                      # aimed at the cluster centre
        L = self.lmp
        L.command(f"create_atoms {TYPES[elem]} single {p[0]:.6f} {p[1]:.6f} {p[2]:.6f} units box")
        new_id = int(np.max(L.numpy.extract_atom("id")[:L.get_natoms()]))
        L.command(f"group newat id {new_id}")
        L.command(f"set group newat charge 0.0")
        L.command(f"velocity newat set {v[0]:.8f} {v[1]:.8f} {v[2]:.8f} units box")
        L.command("group newat delete")
        self.state["recent"][str(new_id)] = self.step
        self.event("insert", elem, int(elem == "C"), int(elem == "O"))

    def remove_desorbed(self):
        ids, types, pos = self.atoms()
        pos, _ = unwrap_about_ni(pos, types, self.box)
        comps, _ = fragments(pos, types)
        is_ni = types == TYPES["Ni"]
        ni = pos[is_ni]
        from scipy.spatial import cKDTree
        tree = cKDTree(ni)
        grace = self.ps(1.0)
        recent = {int(k): s for k, s in self.state["recent"].items() if self.step - s < grace}
        self.state["recent"] = {str(k): s for k, s in recent.items()}
        doomed = []
        for comp in comps:
            if is_ni[comp].any():
                continue
            if any(int(ids[i]) in recent for i in comp):
                continue
            if tree.query(pos[comp])[0].min() > self.cfg["desorb_dist"]:
                doomed.append(comp)
        if not doomed:
            return
        all_ids = []
        for comp in doomed:
            t = types[comp]
            all_ids += [int(ids[i]) for i in comp]
            self.events_pending.append((formula(t), int((t == 1).sum()), int((t == 3).sum())))
        L = self.lmp
        L.command("group gone id " + " ".join(map(str, all_ids)))
        L.command("delete_atoms group gone compress no")
        L.command("group gone delete")
        for f, nC, nO in self.events_pending:
            self.event("desorb", f, nC, nO)
        self.events_pending = []

    def checkpoint(self, final=False):
        L = self.lmp
        tmp = self.f_ckpt + ".tmp"
        L.command(f"write_data {tmp} nocoeff")
        os.replace(tmp, self.f_ckpt)
        if final:
            L.command(f"write_data {os.path.join(self.dir, 'final.data')} nocoeff")
        self.state.update(step=self.step,
                          traj_bytes=os.path.getsize(self.f_traj) if os.path.exists(self.f_traj) else 0,
                          events_bytes=os.path.getsize(self.f_events))
        write_json_atomic(self.f_state, self.state)

    # ------------------------------------------------------------------ main loop
    def run(self):
        self.start_lammps()
        self.events_pending = []
        L = self.lmp
        next_ckpt = (self.step // self.ckpt_every + 1) * self.ckpt_every
        t_last, step_last = time.time(), self.step
        while self.step < self.total:
            j = self.state["j"]
            # ---- events due now
            if j < len(self.seq) and self.step == self.E + j * self.ins_every:
                self.insert(self.seq[j])
                self.state["j"] = j + 1
            if self.cfg["remove_desorbed"] and self.step > 0 and self.step % self.desorb_every == 0:
                self.remove_desorbed()
            if self.step > 0 and self.step % self.recenter_every == 0:
                self.recenter()
            self.state["phase"] = ("equil" if self.step < self.E else
                                   "deposit" if self.step < self.D else "post")
            if self.step >= next_ckpt:
                self.checkpoint()
                next_ckpt += self.ckpt_every
                now = time.time()
                rate = (self.step - step_last) * self.dt / 1000 / max(now - t_last, 1e-9) * 3600
                t_last, step_last = now, self.step
                ids, types, _ = self.atoms()
                self.log(f"t={self.step*self.dt/1000:9.2f} ps  phase={self.state['phase']:7s} "
                         f"inserted={self.state['j']}/{len(self.seq)}  nC={int((types==1).sum())} "
                         f"nO={int((types==3).sum())}  T={L.get_thermo('temp'):7.1f} K  "
                         f"speed={rate:6.1f} ps/h")
                if self.max_sec and time.time() - self.t_wall0 > self.max_sec:
                    self.log("max-hours reached, stopping cleanly at checkpoint")
                    sys.exit(3)
            # ---- run to the next event
            cands = [self.total, next_ckpt,
                     (self.step // self.recenter_every + 1) * self.recenter_every]
            if self.cfg["remove_desorbed"]:
                cands.append((self.step // self.desorb_every + 1) * self.desorb_every)
            j = self.state["j"]
            if j < len(self.seq):
                cands.append(self.E + j * self.ins_every)
            nxt = min(x for x in cands if x > self.step)
            L.command(f"run {nxt - self.step}")
            self.step = nxt
        self.state["phase"] = "done"
        self.checkpoint(final=True)
        self.log(f"finished: t={self.step*self.dt/1000:.2f} ps, wall {(time.time()-self.t_wall0)/3600:.2f} h")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("case_dir")
    ap.add_argument("--threads", type=int, default=int(os.environ.get("OMP_NUM_THREADS", "1")))
    ap.add_argument("--max-hours", type=float, default=None,
                    help="stop cleanly (exit 3) at the first checkpoint after this many hours")
    a = ap.parse_args()
    Case(a.case_dir, a.threads, a.max_hours).run()
