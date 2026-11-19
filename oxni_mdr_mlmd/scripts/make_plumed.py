"""Write PLUMED (OPES-MetaD) inputs for the key elementary steps of a built model.

  python make_plumed.py ../structures/S1_ni211_meo_pore_gas.extxyz --mode ch4 --out ../plumed/S1_ch4.dat
  python make_plumed.py ../structures/S1_ref_ni211_bare_gas.extxyz --mode ch4 --out ../plumed/S1ref_ch4.dat
  python make_plumed.py ../structures/S2_ni211_meo_pore_CO.extxyz --mode coke --out ../plumed/S2_coke.dat

modes
  ch4  : first C-H cleavage (CH4 -> CH3* + H*) of the CH4 nearest the pore axis. CV = C-H
         coordination (~3.8 -> ~2.9). Walls keep the molecule in a cylinder above the pore
         (radius pore-1 A, height `--zmax`): this is the gas reference state, so the rate
         constant refers to one molecule in V_ref and is scaled to the real pressure with
         p*V_ref/kT. CH4 never dissociates in unbiased gas-collision MD (~1e-8 per collision).
  co2  : first C-O cleavage (CO2 -> CO* + O*), CV = C-O coordination of one CO2 (2 -> 1).
  coke : 2-D FES of carbon gasification vs coking among pre-adsorbed C*/O* (S2 models):
         CV1 = N(C-O) (CO formed), CV2 = N(C-C) (C-C bonds = coke nuclei). Lattice-O uptake by
         C (Mars-van Krevelen) near the pore is printed, not biased.
Units: UNITS LENGTH=A ENERGY=eV TIME=fs (PLUMED indices are 1-based). Untested here
(no PLUMED kernel in this container): check with `plumed driver` on a trajectory first.
"""
import argparse
import os
import sys

import numpy as np
from ase.io import read

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
    tags = at.get_tags()
    L = np.diag(at.cell)[:2]
    c = np.array(at.info.get("pore_center_xy", L / 2))
    r = float(at.info.get("pore_radius", 6.0))
    ni = np.where(np.isin(tags, [TAG["ni"], TAG["ni_fixed"]]))[0]
    d = at.positions[ni, :2] - c
    d -= np.round(d / L) * L
    near = ni[np.linalg.norm(d, axis=1) < r]
    z_surf = float(at.positions[near, 2].max())
    return c, r, z_surf


def molecule_near_axis(at, kind, c):
    sym = np.array(at.get_chemical_symbols())
    tags = at.get_tags()
    L = np.diag(at.cell)[:2]
    cs = np.where((sym == "C") & (tags == TAG["gas"]))[0]
    best, out = np.inf, None
    for ci in cs:
        d = at.get_distances(ci, np.where(tags == TAG["gas"])[0], mic=True)
        nb = np.where(tags == TAG["gas"])[0][d < 1.4]
        nb = nb[nb != ci]
        el = set(sym[nb])
        if (kind == "ch4" and el == {"H"} and len(nb) == 4) or (kind == "co2" and el == {"O"} and len(nb) == 2):
            dxy = at.positions[ci, :2] - c
            dxy -= np.round(dxy / L) * L
            score = np.linalg.norm(dxy) + 0.3 * at.positions[ci, 2]
            if score < best:
                best, out = score, (ci, list(nb))
    return out


def header(T):
    return ["UNITS LENGTH=A ENERGY=eV TIME=fs", "# temperature %g K" % T]


def mol_restraints(ci, c, r, z_surf, zmax):
    return [
        "ctr: FIXEDATOM AT=%.3f,%.3f,%.3f" % (c[0], c[1], z_surf),
        "d: DISTANCE ATOMS=ctr,%d COMPONENTS" % (ci + 1),
        "rxy: CUSTOM ARG=d.x,d.y FUNC=sqrt(x*x+y*y) PERIODIC=NO",
        "uwz: UPPER_WALLS ARG=d.z AT=%.2f KAPPA=5.0" % zmax,
        "uwr: UPPER_WALLS ARG=rxy AT=%.2f KAPPA=5.0" % max(r - 1.0, 2.0),
    ]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("structure")
    ap.add_argument("--mode", choices=["ch4", "co2", "coke"], required=True)
    ap.add_argument("--T", type=float, default=873.0)
    ap.add_argument("--barrier", type=float, default=None, help="OPES BARRIER (eV)")
    ap.add_argument("--pace", type=int, default=500)
    ap.add_argument("--zmax", type=float, default=6.5,
                    help="max height of C above exposed Ni (A); with the lateral wall it defines the gas reference volume")
    ap.add_argument("--stride", type=int, default=50)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    at = read(a.structure)
    c, r, z_surf = pore_geometry(at)
    L = header(a.T)
    if a.mode in ("ch4", "co2"):
        ci, nb = molecule_near_axis(at, a.mode, c)
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
        L += mol_restraints(ci, c, r, z_surf, a.zmax)
        v_ref = np.pi * max(r - 1.0, 2.0) ** 2 * a.zmax
        L.insert(2, "# gas reference volume V_ref ~ %.0f A^3: events per site per s = k(V_ref) * p*V_ref/kT / N_exposed_Ni" % v_ref)
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
    with open(a.out, "w") as fh:
        fh.write("# generated by make_plumed.py from %s (mode %s)\n" % (os.path.basename(a.structure), a.mode))
        fh.write("\n".join(L) + "\n")
    print(open(a.out).read())


if __name__ == "__main__":
    main()
