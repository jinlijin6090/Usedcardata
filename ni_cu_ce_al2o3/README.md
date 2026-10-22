# NiCu / Ce-γ-Al2O3 model catalyst

This folder builds and renders an atomistic model of a Ni particle on γ-Al2O3. Cu replaces
the Ni atoms on the low-coordination (corner and edge) sites, Ce replaces some surface Al,
and some surface O is removed to make O vacancies. It also has three analyses:

1. The Cu/Ni ratio needed to put Cu on every corner and edge site, as a function of particle size.
2. How Ce doping of Al2O3 changes the O-vacancy energetics.
3. How particle size changes CH4-activation activity and coking.

## Layout

| path | content |
|---|---|
| `scripts/model_lib.py` | builders: defective-spinel γ-Al2O3(100) slab, truncated-octahedral Ni particles, site classification, Cu decoration, Ce doping + V_O |
| `scripts/build_and_render.py` | builds the models, writes `structures/*.xyz/.cif`, renders `renders/*.png` with POV-Ray |
| `scripts/cu_segregation.py` | EMT segregation energies of Cu in Ni + site statistics for 1.6–9.8 nm particles → `data/cu_segregation.json` |
| `scripts/ce_vacancy_mace.py` | MACE-MP-0 (medium) O-vacancy formation energies on Ce-doped γ-Al2O3(100) → `data/ce_vacancy_mace.json` |
| `scripts/fig_*.py` | figures 1–4 in `figures/` |

Run order: `build_and_render.py`, `cu_segregation.py`, `ce_vacancy_mace.py` (about 1 h on 4 CPU cores),
then `fig_renders.py`, `fig_cu.py`, `fig_ce.py`, `fig_activity_coke.py`.
Requirements: ase, numpy, scipy, matplotlib, povray, torch, mace-torch.

## Model details and assumptions

**γ-Al2O3.** The MgAl2O4 spinel lattice (a = 7.911 Å) has all its cation sites filled with Al.
Cation vacancies are then added at random (mostly octahedral sites, sub-surface and bulk)
until Al:O = 2:3. The (100) surface ends in the mixed octahedral-Al/O layer. The main slab is
5×5 cells (39.6 Å × 39.6 Å) with 5 mixed layers, Al600O900 before doping.

**Ni particle.** The particle is an fcc truncated octahedron (a = 3.524 Å) sitting on its
(100) facet. Bottom layers are cut off to mimic wetting. The bottom layer counts as the
metal/oxide interface. Exposed sites are classified by metal coordination number:
corner CN ≤ 6, edge CN = 7, (100) CN = 8, (111) CN = 9.

**Cu.** In the renders every corner and edge atom is Cu (the T = 0 limit). The quantitative
model uses EMT segregation energies of one Cu atom in a 586-atom Ni particle:
corner −0.52, edge −0.49, (100) −0.51, (111) −0.36, sub-surface −0.18 eV. These go into a
Fermi–Dirac (McLean) site-occupancy model at 973 K, with the Cu chemical potential solved
from the total Cu/Ni ratio. Cu–Cu interactions are ignored. A bond-cutting variant
(segregation energy linear in missing bonds) is plotted as a sensitivity check.

**Ce and V_O.** Ce replaces surface octahedral Al. The isovalent Ce_Al substitution keeps the
cell neutral. O vacancies are removed next to Ce. Formation energies are
E_vf = E(slab−O) + ½E(O2) − E(slab), computed with the universal ML potential MACE-MP-0.
MACE-MP-0 is trained on PBE/PBE+U data with no U on Ce, so Ce-4f localisation (Ce4+→Ce3+) is
not described properly. Read the trend relative to undoped Al2O3, not the absolute numbers.
The equilibrium vacancy enhancement assumes random Ce on the surface: each surface O has
two surface-cation neighbours, the Ce-neighbour count is binomial, and occupancy is
Boltzmann-weighted.

**Activity and coking (descriptor model, arbitrary units).**
- CH4 dissociation barriers: corner 0.75, edge 0.80, (100) 0.90, (111) 1.00, perimeter 0.85 eV,
  and 1.45 eV on sites covered by Cu. These are representative values from the Ni DFT literature.
- Carbon nucleation is proportional to the rate on free Ni step/corner sites (Bengaard et al.,
  J. Catal. 2002), multiplied by a logistic size factor (d_c = 4 nm, width 0.8 nm) for the
  suppression of graphene/whisker growth on small particles.
- Carbon gasification scales with the number of perimeter sites × √(V_O enhancement from Ce).
- Coke index = nucleation / (nucleation + gasification).

These are qualitative structure–property models for building hypotheses and making
illustrations. Before you quote numbers, validate them with DFT(+U) and experiment
(H2-TPR, TGA/TPO, TEM size distributions).
