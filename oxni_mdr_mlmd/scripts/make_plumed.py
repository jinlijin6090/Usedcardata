"""Write PLUMED (OPES-MetaD) inputs for elementary steps of a built model.

  python make_plumed.py ../structures/S1_ni211_meo_pore_nogas.extxyz --mode ch4 --out ../plumed/S1_ch4_873K.dat
  python make_plumed.py ../structures/S1_ref_ni211_bare_nogas.extxyz --mode ch4 --out ../plumed/S1ref_ch4_873K.dat
  python make_plumed.py ../structures/S2_ni211_meo_pore_CO.extxyz --mode coke --out ../plumed/S2_coke_873K.dat

modes
  ch4  : first C-H cleavage (CH4 -> CH3* + H*) of ONE CH4. CV = C-H coordination (~3.8 -> ~2.9).
  co2  : first C-O cleavage (CO2 -> CO* + O*) of ONE CO2. CV = C-O coordination (~1.9 -> ~1.0).
         For both, walls keep the molecule in a reference cylinder on the pore axis: radius
         pore-1 A, from the exposed Ni up to `--zmax`. Default height: 8 A above a bare
         surface; (film top) + 6 A when an oxide film is present, so that covered and bare
         models share the same reference state, gas just outside the solid. The script puts
         the molecule into the cylinder (it builds one if the structure has none, or moves the
         nearest one) and writes that structure next to the input as *_start.extxyz: run the
         MD from that file. Use gas-free models; other gas molecules would crowd the cylinder.
         CH4 never dissociates in unbiased gas-collision MD (~1e-8 per collision at 873 K).
  coke : 2-D FES among pre-adsorbed C*/O* (S2 models): CV1 = N(C-O), CV2 = N(C-C). Lattice-O
         uptake by C near the pore is printed, not biased. MECHANISM_PLAN.md no longer uses
         this as a core test (C2 formation is not coke); kept as a template.

Rate constants: the run gives k for one molecule confined to the cylinder. Per-site rate at
pressure p = k * p * V_eff / kT / N_exposed_Ni with V_eff = integral over the cylinder of
exp(-w/kT), w = potential of mean force of the intact molecule relative to the flat gas part
(reweighted d.z histogram). V_eff equals the geometric volume only for a weakly bound molecule.

Not available here: inputs for C*+O*, CH*+O*, CO dissociation (direct / via COH), OPES-flooding
and committor runs needed by work package A3 of MECHANISM_PLAN.md.
Units: UNITS LENGTH=A ENERGY=eV TIME=fs (PLUMED indices are 1-based). Never run with PLUMED
(no PLUMED kernel in this container): check with `plumed driver` on a trajectory first.
"""
import argparse
import os
import sys

import numpy as np
from ase.build import molecule
from ase.io import read, write

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from meo_lib import TAG  # noqa: E402


def sw(r0, nn=6, mm=12):
    """Rational switching function s(r) = (1-(r/r0)^nn)/(1-(r/r0)^mm), s(r0) = 0.5.

    r0 is put near the transition-state bond length, so reactant, transition state and
    product are all resolved and the bias still acts at neighbour-site distances (2.5-3 A).
    (A narrow D_0/R_0 switch squeezes the transition state onto the product value.)
    """
    return "SWITCH={RATIONAL R_0=%.2f NN=%d MM=%d}" % (r0, nn, mm)


def lst(idx):
    return ",".join(str(i + 1) for i in idx)


def pore_geometry(at):
    """Pore axis (x, y), pore radius, height of the exposed Ni, top of the oxide film (or None)."""
    tags = at.get_tags()
    L = np.diag(at.cell)[:2]
    c = np.array(at.info.get("pore_center_xy", L / 2))
    r = float(at.info.get("pore_radius", 6.0))
    ni = np.where(np.isin(tags, [TAG["ni"], TAG["ni_fixed"]]))[0]
    d = at.positions[ni, :2] - c
    d -= np.round(d / L) * L
    near = ni[np.linalg.norm(d, axis=1) < r]
    z_surf = float(at.positions[near, 2].max())
    film = np.where(tags == TAG["meo"])[0]
    z_film = float(at.positions[film, 2].max()) if len(film) else None
    return c, r, z_surf, z_film


def molecule_near_axis(at, kind, c):
    """Index of C and of its ligands for the intact gas molecule closest to the pore axis."""
    sym = np.array(at.get_chemical_symbols())
    tags = at.get_tags()
    L = np.diag(at.cell)[:2]
    gas = np.where(tags == TAG["gas"])[0]
    best, out = np.inf, None
    for ci in gas[sym[gas] == "C"]:
        d = at.get_distances(ci, gas, mic=True)
        nb = gas[d < 1.4]
        nb = nb[nb != ci]
        el = set(sym[nb])
        if (kind == "ch4" and el == {"H"} and len(nb) == 4) or (kind == "co2" and el == {"O"} and len(nb) == 2):
            dxy = at.positions[ci, :2] - c
            dxy -= np.round(dxy / L) * L
            score = np.linalg.norm(dxy) + 0.3 * at.positions[ci, 2]
            if score < best:
                best, out = score, (int(ci), [int(i) for i in nb])
    return out


def in_cylinder(at, ci, c, z_surf, r_wall, zmax):
    L = np.diag(at.cell)[:2]
    dxy = at.positions[ci, :2] - c
    dxy -= np.round(dxy / L) * L
    dz = at.positions[ci, 2] - z_surf
    return np.linalg.norm(dxy) < r_wall - 0.3 and 1.5 < dz < zmax - 0.3


def put_molecule_in_cylinder(at, kind, c, z_surf, z_film, r_wall, zmax, min_dist=2.6):
    """Return (atoms, ci, nb, changed): one intact molecule with its C on the pore axis.

    Start height: in the gas part of the cylinder (above the film if there is one). The
    molecule nearest the axis is moved rigidly; if the structure has none, one is built.
    The cell is made tall enough that the cylinder top is >= 8 A below the periodic image,
    and info['wall_z'] (reflecting wall of run_md.py) is lifted above the PLUMED wall.
    """
    at = at.copy()
    changed = False
    found = molecule_near_axis(at, kind, c)
    if found and in_cylinder(at, found[0], c, z_surf, r_wall, zmax):
        ci, nb = found
    else:
        if found:
            ci, nb = found
            ids = [ci] + nb
            rel = at.get_distances(ci, ids, mic=True, vector=True)  # whole molecule around C
        else:
            m = molecule("CH4" if kind == "ch4" else "CO2")
            rel = m.positions - m.positions[m.get_chemical_symbols().index("C")]
            order = [m.get_chemical_symbols().index("C")] + [i for i, s in enumerate(m.get_chemical_symbols()) if s != "C"]
            rel = rel[order]
            m = m[order]
            m.set_tags(np.full(len(m), TAG["gas"]))
            ci = len(at)
            nb = list(range(ci + 1, ci + len(m)))
            ids = [ci] + nb
            at += m
        others = np.setdiff1d(np.arange(len(at)), ids)
        h0 = (z_film - z_surf) + 3.5 if z_film is not None else 5.5   # preferred start height
        trial = sorted(np.arange(h0 - 1.0, zmax - 0.8, 0.25), key=lambda h: abs(h - h0))
        for h in trial:
            pos = np.r_[c, z_surf + h] + rel
            tmp = at.copy()
            tmp.positions[ids] = pos
            dmin = min(tmp.get_distances(i, others, mic=True).min() for i in ids) if len(others) else 99.0
            if dmin > min_dist:
                at.positions[ids] = pos
                break
        else:
            raise RuntimeError("no free start position in the reference cylinder; raise --zmax")
        changed = True
    need = z_surf + zmax + 8.0
    if at.cell[2, 2] < need:
        at.cell[2, 2] = need
        changed = True
    wall = z_surf + zmax + 3.0   # run_md.py's reflecting wall must stay above the PLUMED wall
    if float(at.info.get("wall_z", wall)) < wall:
        at.info["wall_z"] = float(wall)
        changed = True
    return at, ci, nb, changed


def header(T):
    return ["UNITS LENGTH=A ENERGY=eV TIME=fs", "# temperature %g K" % T]


def mol_restraints(ci, c, r_wall, z_surf, zmax):
    return [
        "ctr: FIXEDATOM AT=%.3f,%.3f,%.3f" % (c[0], c[1], z_surf),
        "d: DISTANCE ATOMS=ctr,%d COMPONENTS" % (ci + 1),
        "rxy: CUSTOM ARG=d.x,d.y FUNC=sqrt(x*x+y*y) PERIODIC=NO",
        "uwz: UPPER_WALLS ARG=d.z AT=%.2f KAPPA=5.0" % zmax,
        "uwr: UPPER_WALLS ARG=rxy AT=%.2f KAPPA=5.0" % r_wall,
    ]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("structure")
    ap.add_argument("--mode", choices=["ch4", "co2", "coke"], required=True)
    ap.add_argument("--T", type=float, default=873.0)
    ap.add_argument("--barrier", type=float, default=None, help="OPES BARRIER (eV)")
    ap.add_argument("--pace", type=int, default=500)
    ap.add_argument("--zmax", type=float, default=None,
                    help="top of the reference cylinder above the exposed Ni (A); "
                         "default 8 (bare) or film top + 6 (with an oxide film)")
    ap.add_argument("--stride", type=int, default=50)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    at = read(a.structure)
    c, r, z_surf, z_film = pore_geometry(at)
    L = header(a.T)
    start = None
    if a.mode in ("ch4", "co2"):
        r_wall = max(r - 1.0, 2.0)
        zmax = a.zmax if a.zmax is not None else ((z_film - z_surf) + 6.0 if z_film is not None else 8.0)
        at, ci, nb, changed = put_molecule_in_cylinder(at, a.mode, c, z_surf, z_film, r_wall, zmax)
        if changed:
            start = os.path.splitext(os.path.abspath(a.out))[0] + "_start.extxyz"
        if a.mode == "ch4":
            # CH4 ~3.8, transition state (one C-H ~1.6 A) ~3.4, CH3* + H* ~2.9
            L += ["cn: COORDINATION GROUPA=%d GROUPB=%s %s" % (ci + 1, lst(nb), sw(1.60, 8, 16)),
                  "lw: LOWER_WALLS ARG=cn AT=2.5 KAPPA=20.0   # stop at CH3*+H* (first C-H cleavage)"]
            barrier = a.barrier or 1.5
        else:
            # CO2 ~1.9, transition state (one C-O ~1.85 A) ~1.45, CO* + O* ~1.0
            L += ["cn: COORDINATION GROUPA=%d GROUPB=%s %s" % (ci + 1, lst(nb), sw(1.85, 8, 16)),
                  "lw: LOWER_WALLS ARG=cn AT=0.6 KAPPA=20.0   # stop at CO*+O*"]
            barrier = a.barrier or 1.5
        L += mol_restraints(ci, c, r_wall, z_surf, zmax)
        v_geo = np.pi * r_wall ** 2 * zmax
        L.insert(2, "# reference cylinder on the pore axis: radius %.1f A, 0-%.1f A above the exposed Ni%s; geometric volume %.0f A^3"
                 % (r_wall, zmax, " (film top at %.1f A)" % (z_film - z_surf) if z_film is not None else "", v_geo))
        L.insert(3, "# per-site rate at pressure p = k * p*V_eff/kT / N_exposed_Ni, V_eff = integral of exp(-w(z)/kT) over the cylinder (reweighted d.z histogram)")
        L += ["opes: OPES_METAD ARG=cn PACE=%d BARRIER=%.2f TEMP=%g FILE=KERNELS STATE_WFILE=STATE STATE_WSTRIDE=10000"
              % (a.pace, barrier, a.T),
              "PRINT ARG=cn,d.z,rxy,opes.bias,uwz.bias,uwr.bias,lw.bias STRIDE=%d FILE=COLVAR" % a.stride]
    else:
        sym = np.array(at.get_chemical_symbols())
        tags = at.get_tags()
        C = np.where((sym == "C") & (tags == TAG["ads"]))[0]
        O = np.where((sym == "O") & (tags == TAG["ads"]))[0]
        Lxy = np.diag(at.cell)[:2]
        dxy = at.positions[:, :2] - c
        dxy -= np.round(dxy / Lxy) * Lxy
        Olat = np.where((sym == "O") & (tags == TAG["meo"]) & (np.linalg.norm(dxy, axis=1) < r + 6.0))[0]
        # bonded pair ~0.9, transition state ~0.5, neighbouring sites (2.5 A) ~0.2, far apart ~0
        L += ["nco: COORDINATION GROUPA=%s GROUPB=%s %s" % (lst(C), lst(O), sw(1.90)),
              "ncc: COORDINATION GROUPA=%s %s" % (lst(C), sw(2.00))]
        if len(Olat):
            L += ["ncolat: COORDINATION GROUPA=%s GROUPB=%s %s" % (lst(C), lst(Olat), sw(1.90))]
        barrier = a.barrier or 2.0
        L += ["opes: OPES_METAD ARG=nco,ncc PACE=%d BARRIER=%.2f TEMP=%g FILE=KERNELS STATE_WFILE=STATE STATE_WSTRIDE=10000"
              % (a.pace, barrier, a.T),
              "PRINT ARG=nco,ncc%s,opes.bias STRIDE=%d FILE=COLVAR" % (",ncolat" if len(Olat) else "", a.stride)]
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    src = os.path.basename(start) if start else os.path.basename(a.structure)
    with open(a.out, "w") as fh:
        fh.write("# generated by make_plumed.py (mode %s); run the MD from %s\n" % (a.mode, src))
        fh.write("\n".join(L) + "\n")
    if start:
        write(start, at)
        print("start structure with the molecule inside the reference cylinder:", start)
    print(open(a.out).read())


if __name__ == "__main__":
    main()
