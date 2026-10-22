"""Build the NiCu/Ce-gamma-Al2O3 models, save structures, and render them with POV-Ray."""
import os
import numpy as np
from ase import Atoms
from ase.io import write

from model_lib import (COLORS, RADII, classify_sites, decorate_cu, dope_ce_and_vacancies,
                       gamma_alumina_slab, ni_particle, place_on)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
STRUCT = os.path.join(ROOT, "structures")
RENDER = os.path.join(ROOT, "renders")


def render(atoms, name, rotation, width=1800, vac_pos=None, transparent_support=False,
           highlight=None):
    a = atoms.copy()
    if vac_pos is not None and len(vac_pos):
        a += Atoms("X" * len(vac_pos), positions=vac_pos)
    sym = a.get_chemical_symbols()
    colors = [COLORS[s] for s in sym]
    radii = [RADII[s] for s in sym]
    trans = [0.35 if s == "X" else 0.0 for s in sym]
    if transparent_support:
        trans = [0.65 if s in ("Al", "O") else t for s, t in zip(sym, trans)]
    if highlight is not None:  # dim everything except highlighted indices
        trans = [t if i in highlight else max(t, 0.75) for i, t in enumerate(trans)]
    textures = ["vmd" if s in ("Ni", "Cu", "Ce") else "jmol" for s in sym]
    textures = ["glass" if s == "X" else t for s, t in zip(sym, textures)]
    a.cell = np.zeros((3, 3))  # no unit-cell box in the picture
    pov = os.path.join(RENDER, name + ".pov")
    os.chdir(RENDER)  # POV-Ray resolves the .pov named in the .ini relative to cwd
    write(pov, a, rotation=rotation, radii=radii, colors=colors,
          povray_settings=dict(canvas_width=width, transmittances=trans, textures=textures,
                               celllinewidth=0.0, background="White", transparent=False,
                               area_light=[(2., 3., 40.), "White", .7, .7, 4, 4],
                               point_lights=[((-40., 25., 60.), (0.45, 0.45, 0.45)), ((40., -10., 60.), (0.25, 0.25, 0.25))])
          ).render()
    for ext in (".pov", ".ini"):
        os.remove(pov.replace(".pov", ext))
    return os.path.join(RENDER, name + ".png")


def main_model():
    slab = gamma_alumina_slab(5, 5, n_mixed=5, seed=0)
    part = ni_particle(9, 3, cut_bottom_layers=1)  # ~1.8 nm truncated octahedron
    cn, types = classify_sites(part)
    part = decorate_cu(part, types, fill=("corner", "edge"))
    centre = slab.cell[0, :2] / 2 + slab.cell[1, :2] / 2
    slab, vac, n_ce = dope_ce_and_vacancies(slab, n_ce=14, vac_per_ce=0.6, center_xy=centre,
                                            exclude_r=9.5, seed=3)
    model = place_on(slab, part)
    model.set_array("site", np.r_[np.full(len(slab), "support", object), types].astype(str))
    write(os.path.join(STRUCT, "NiCu_Ce-gAl2O3_main.xyz"), model)
    write(os.path.join(STRUCT, "NiCu_Ce-gAl2O3_main.cif"), model)
    np.savetxt(os.path.join(STRUCT, "NiCu_Ce-gAl2O3_main_Ovacancy_positions.txt"), vac,
               header="x y z (A) of removed surface O (oxygen vacancies)")
    print("main:", model.get_chemical_formula(), "Ce:", n_ce, "V_O:", len(vac),
          "Cu:", (model.symbols == "Cu").sum(), "Ni:", (model.symbols == "Ni").sum())

    render(model, "01_overview_oblique", "35z,-58x", vac_pos=vac)
    render(model, "02_top_view", "0x,0y,0z", vac_pos=vac)
    # close-up of the particle only (Cu on corners/edges)
    p = model[len(slab):]
    render(p, "03_particle_closeup", "35z,-60x", width=1400)
    return model, slab, vac


def size_series():
    """Same Cu-on-edges motif on three particle sizes (same support, same scale)."""
    out = []
    specs = [("small", 5, 1, 0), ("medium", 9, 3, 1), ("large", 13, 4, 2)]
    for tag, L, c, cb in specs:
        slab = gamma_alumina_slab(6, 6, n_mixed=3, seed=0)
        part = ni_particle(L, c, cut_bottom_layers=cb)
        cn, types = classify_sites(part)
        part = decorate_cu(part, types)
        centre = slab.cell[0, :2] / 2 + slab.cell[1, :2] / 2
        slab, vac, _ = dope_ce_and_vacancies(slab, 16, 0.6, centre, exclude_r=16, seed=4)
        m = place_on(slab, part)
        ncu = (part.symbols == "Cu").sum()
        nni = (part.symbols == "Ni").sum()
        d = np.ptp(part.positions[:, 0]) / 10 + 0.25
        print(f"{tag}: {len(part)} atoms, d~{d:.1f} nm, Cu={ncu}, Ni={nni}, Cu/Ni={ncu/nni:.3f}")
        write(os.path.join(STRUCT, f"size_{tag}.xyz"), m)
        out.append(render(m, f"04_size_{tag}", "35z,-58x", width=1200, vac_pos=vac))
    return out


if __name__ == "__main__":
    os.makedirs(STRUCT, exist_ok=True)
    os.makedirs(RENDER, exist_ok=True)
    main_model()
    size_series()
