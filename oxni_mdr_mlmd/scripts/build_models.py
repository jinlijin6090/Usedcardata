"""Build the OxNi / reference model set -> structures/*.extxyz and renders/*.png.

  S0_*  small smoke-test slab (CPU-sized)
  S1_*  Ni(211) + 1 nm MEO film with a 1.2 nm pore + CH4/CO2 gas    | bare Ni(211) reference
  S2_*  same slabs, gas-free, C*/O* pre-loaded in a 1.2 nm disc      | (coke vs gasification)
  S3_*  Ni(111) terrace version of S1
  P1_*  3 nm Ni particle + porous 1 nm MEO shell on Mg(Al)O(100)     | bare particle reference
  P2_*  two 2.5 nm particles side by side (sintering / coalescence)  | bare pair reference

Usage: python build_models.py [--only S0,S1] [--no-render]
"""
import argparse
import os
import sys

import numpy as np
from ase.io import write

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import meo_lib as ml  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SDIR, RDIR = os.path.join(ROOT, "structures"), os.path.join(ROOT, "renders")


def slab_pair(facet, nx, ny, layers, px, py, film_layers, pore_r, seed, o_deficit=0.0):
    ni = ml.ni_slab(facet, nx, ny, layers)
    film = ml.rocksalt_film(px, py, film_layers, ml.MEO_EQUIMOLAR, o_deficit=o_deficit, seed=seed)
    covered, strain = ml.stack_film(ni, film)
    covered = ml.carve_pore(covered, pore_r)
    covered.info["film_strain_xy"] = tuple(map(float, strain))
    ni.info["pore_center_xy"] = covered.info["pore_center_xy"]
    ni.info["pore_radius"] = covered.info["pore_radius"]
    return covered, ni


def build(only=None):
    models = {}
    want = (lambda k: True) if not only else (lambda k: any(k.startswith(o) for o in only))

    if want("S0"):
        cov, bare = slab_pair("211", 9, 8, 4, 6, 7, 5, 5.0, seed=1)
        models["S0_smoke_ni211_meo_pore_gas"] = ml.add_gas(cov, 6, 6, height=14.0, seed=1)
        models["S0_smoke_ref_ni211_bare_gas"] = ml.add_gas(bare, 6, 6, height=14.0, seed=1)

    if want("S1") or want("S2"):
        # Ni(211) 24.4 x 27.4 A, 5 step units (~10 A); film 8 x 9 cells, 6 layers (~1.05 nm),
        # strain +2.8 / +2.6 %; pore radius 6 A exposes ~1.2 nm of stepped Ni
        cov, bare = slab_pair("211", 12, 11, 5, 8, 9, 6, 6.0, seed=2)
        if want("S1"):
            # gas-free versions: anneal the MEO first (protocols/meo_anneal.json), then add
            # gas / adatoms to the annealed frame with add_species.py
            for name, base in (("S1_ni211_meo_pore_nogas", cov), ("S1_ref_ni211_bare_nogas", bare)):
                m = base.copy()
                m.info["wall_z"] = float(m.positions[:, 2].max() + 8.0)
                models[name] = m
            models["S1_ni211_meo_pore_gas"] = ml.add_gas(cov, 12, 12, height=24.0, seed=2)
            models["S1_ref_ni211_bare_gas"] = ml.add_gas(bare, 12, 12, height=24.0, seed=2)
        if want("S2"):
            c, r = cov.info["pore_center_xy"], 5.0
            load = {"C": 4, "O": 4}
            for name, base in (("S2_ni211_meo_pore_CO", cov), ("S2_ref_ni211_bare_CO", bare)):
                m = ml.add_adatoms(base, load, c, r, seed=3)
                m.info["wall_z"] = float(m.positions[:, 2].max() + 10.0)  # keeps desorbed CO near
                m.cell[2, 2] = max(m.cell[2, 2], m.info["wall_z"] + 4.0)
                models[name] = m

    if want("S3"):
        cov, bare = slab_pair("111", 11, 12, 5, 9, 9, 6, 6.0, seed=4)
        models["S3_ni111_meo_pore_gas"] = ml.add_gas(cov, 12, 12, height=24.0, seed=4)
        models["S3_ref_ni111_bare_gas"] = ml.add_gas(bare, 12, 12, height=24.0, seed=4)

    if want("P1"):
        spec = dict(xy=(0.5, 0.5), d=30.0, shell=10.0, n_pores=6, pore_radius=4.5)
        models["P1_np3nm_meo_on_MgAlO"] = ml.particles_on_support([spec], (24, 24), seed=5)
        models["P1_ref_np3nm_bare_on_MgAlO"] = ml.particles_on_support(
            [dict(spec, shell=0)], (24, 24), seed=5)

    if want("P2"):
        px, py = 40, 20
        lx = px * ml.A_MGO / np.sqrt(2)
        d, shell, sgap = 25.0, 10.0, 2.0
        # shells (or Ni surfaces for the bare pair) start ~3 A apart
        half = (d / 2 + sgap + shell + 1.5) / lx
        half_bare = (d / 2 + 1.5) / lx
        pair = [dict(xy=(0.5 - half, 0.5), d=d, shell=shell, n_pores=5),
                dict(xy=(0.5 + half, 0.5), d=d, shell=shell, n_pores=5)]
        pair_bare = [dict(xy=(0.5 - half_bare, 0.5), d=d, shell=0),
                     dict(xy=(0.5 + half_bare, 0.5), d=d, shell=0)]
        models["P2_pair_np2.5nm_meo"] = ml.particles_on_support(pair, (px, py), seed=6)
        models["P2_ref_pair_np2.5nm_bare"] = ml.particles_on_support(pair_bare, (px, py), seed=6)
    return models


def render(name, atoms):
    col = ml.colors_for(atoms)
    kw = dict(colors=col, radii=0.85, scale=12 if len(atoms) < 3000 else 6)
    if name.startswith("P"):
        c = atoms.positions[atoms.get_tags() == ml.TAG["ni"]].mean(axis=0)
        cut = atoms[np.abs(atoms.positions[:, 1] - c[1]) < 4.0]  # cross-section, like Fig. 1i
        write(os.path.join(RDIR, name + "_section.png"), cut, rotation="-90x",
              **dict(kw, colors=ml.colors_for(cut)))
        write(os.path.join(RDIR, name + "_top.png"), atoms, **kw)
    else:
        write(os.path.join(RDIR, name + "_side.png"), atoms, rotation="-80x,10y", **kw)
        top = atoms[atoms.get_tags() != ml.TAG["gas"]]
        top.cell = atoms.cell
        write(os.path.join(RDIR, name + "_top.png"), top, **dict(kw, colors=ml.colors_for(top)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default="")
    ap.add_argument("--no-render", action="store_true")
    a = ap.parse_args()
    os.makedirs(SDIR, exist_ok=True)
    os.makedirs(RDIR, exist_ok=True)
    models = build([o for o in a.only.split(",") if o])
    with open(os.path.join(SDIR, "SUMMARY.txt"), "a" if a.only else "w") as fh:
        for name, at in models.items():
            info = {k: v for k, v in at.info.items() if k != "adsorbate_info"}
            at.info = info
            write(os.path.join(SDIR, name + ".extxyz"), at)
            txt = "%s\n%s\n  info: %s\n" % (name, ml.summary(at), info)
            print(txt)
            fh.write(txt + "\n")
            if not a.no_render:
                render(name, at)


if __name__ == "__main__":
    main()
