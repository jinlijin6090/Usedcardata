"""Builders for Cu-decorated Ni particles on Ce-doped gamma-Al2O3(100) with O vacancies.

gamma-Al2O3 is modelled as a defective spinel (Fd-3m, a = 7.911 A): the Mg (8a, tetrahedral)
and Al (16d, octahedral) sites of MgAl2O4 are all filled with Al, then cation vacancies are
introduced until the Al:O ratio is 2:3 (21 1/3 Al per 32 O). The slab is terminated by a
mixed octahedral-Al/O (100) layer, which is the common Digne-type termination.
"""
import numpy as np
from ase import Atoms
from ase.build import surface
from ase.neighborlist import neighbor_list
from ase.spacegroup import crystal
from scipy.spatial import cKDTree

A_SPINEL = 7.911
A_NI = 3.524

# colours (RGB 0-1) and radii (A) used by every render
COLORS = {
    "O": (0.86, 0.16, 0.14),
    "Al": (0.80, 0.78, 0.86),
    "Ce": (1.00, 0.78, 0.10),
    "Ni": (0.30, 0.72, 0.42),
    "Cu": (0.85, 0.45, 0.18),
    "X": (0.10, 0.85, 1.00),  # oxygen vacancy marker
}
RADII = {"O": 0.88, "Al": 0.66, "Ce": 1.08, "Ni": 1.22, "Cu": 1.26, "X": 0.95}


def spinel_bulk():
    b = crystal(["Mg", "Al", "O"], [(1 / 8,) * 3, (1 / 2,) * 3, (0.2625,) * 3],
                spacegroup=227, setting=2, cellpar=[A_SPINEL] * 3 + [90] * 3)
    site = np.where(np.array(b.get_chemical_symbols()) == "Mg", 1, 0)  # 1 = tetrahedral
    b.set_array("tet", site)
    b.symbols[b.symbols == "Mg"] = "Al"
    return b


def gamma_alumina_slab(nx, ny, n_mixed=5, vacuum=12.0, seed=0):
    """(100) slab with `n_mixed` octahedral-Al/O layers (tetrahedral layers in between)."""
    rng = np.random.default_rng(seed)
    s = surface(spinel_bulk(), (1, 0, 0), 2, vacuum=None, periodic=True)
    s = s.repeat((nx, ny, 1))
    z = s.positions[:, 2]
    is_o = s.symbols == "O"
    # mixed layers are where O sits; keep n_mixed of them plus the tetrahedral layers between
    o_levels = np.unique(np.round(z[is_o], 1))
    mixed = []
    for zl in o_levels:  # O planes come in pairs ~0.2 A apart: merge them
        if not mixed or zl - mixed[-1][-1] > 0.5:
            mixed.append([zl])
        else:
            mixed[-1].append(zl)
    centers = [np.mean(m) for m in mixed][:n_mixed]
    zmin, zmax = centers[0] - 0.5, centers[-1] + 0.5
    s = s[(z > zmin) & (z < zmax)]
    s.positions[:, 2] -= s.positions[:, 2].min()
    s.cell[2, 2] = s.positions[:, 2].max() + vacuum
    s.pbc = (True, True, False)

    # cation vacancies -> Al2O3 stoichiometry. Keep the top (surface) layer intact,
    # remove mostly octahedral Al from the sub-surface/bulk layers (~1:3 tet:oct).
    n_o = (s.symbols == "O").sum()
    n_al = (s.symbols == "Al").sum()
    n_remove = n_al - int(round(2 * n_o / 3))
    z = s.positions[:, 2]
    al = np.where((s.symbols == "Al") & (z < z.max() - 1.5) & (z > 1.5))[0]
    tet = s.get_array("tet")
    w = np.where(tet[al] == 1, 1.0, 1.5)
    drop = rng.choice(al, size=n_remove, replace=False, p=w / w.sum())
    del s[drop]
    return s


def surface_layers(slab):
    z = slab.positions[:, 2]
    top = z.max()
    return top, (z > top - 0.6)


def dope_ce_and_vacancies(slab, n_ce, vac_per_ce=0.5, center_xy=None, exclude_r=0.0,
                          near_r=None, seed=1):
    """Replace surface octahedral Al by Ce and remove O neighbours of Ce (O vacancies).

    Ce4+ on an Al3+ site that is reduced to Ce3+ is what makes the neighbouring O easy to
    remove, so vacancies are placed next to Ce. Returns (slab, vacancy_positions).
    """
    rng = np.random.default_rng(seed)
    slab = slab.copy()
    top, is_top = surface_layers(slab)
    xy = slab.positions[:, :2]
    cand = np.where(is_top & (slab.symbols == "Al"))[0]
    if center_xy is not None:
        d = np.linalg.norm(xy[cand] - center_xy, axis=1)
        keep = d > exclude_r
        if near_r is not None:
            keep &= d < near_r
        cand = cand[keep]
    # spread Ce out: greedy pick with minimum separation
    rng.shuffle(cand)
    ce = []
    for i in cand:
        if all(slab.get_distance(i, j, mic=True) > 5.0 for j in ce):
            ce.append(i)
        if len(ce) == n_ce:
            break
    slab.symbols[ce] = "Ce"
    # Ce is bigger than Al: push it slightly outwards
    slab.positions[ce, 2] += 0.35

    n_vac = int(round(vac_per_ce * len(ce)))
    i_idx, j_idx = neighbor_list("ij", slab, 2.6)
    vac = []
    for c in rng.permutation(ce):
        if len(vac) >= n_vac:
            break
        o_nb = [j for i, j in zip(i_idx, j_idx) if i == c and slab.symbols[j] == "O"
                and slab.positions[j, 2] > top - 1.0 and j not in vac]
        if o_nb:
            vac.append(int(rng.choice(o_nb)))
    vac_pos = slab.positions[vac].copy()
    del slab[vac]
    return slab, vac_pos, len(ce)


def truncated_octahedron(length, cutoff, a=A_NI):
    """fcc truncated octahedron, same convention as ase.cluster.Octahedron(length, cutoff).

    Integer lattice in units of a/2: fcc points have i+j+k even (or odd, chosen so that
    the octahedron vertex is a lattice point). {111} faces: |i|+|j|+|k| <= L-1,
    {100} faces: max(|i|,|j|,|k|) <= L-1-cutoff.
    """
    n111 = length - 1
    n100 = length - 1 - cutoff
    r = np.arange(-n100, n100 + 1)
    i, j, k = np.meshgrid(r, r, r, indexing="ij")
    i, j, k = i.ravel(), j.ravel(), k.ravel()
    keep = ((i + j + k) % 2 == n111 % 2) & (np.abs(i) + np.abs(j) + np.abs(k) <= n111)
    pos = np.c_[i[keep], j[keep], k[keep]] * (a / 2)
    return Atoms("Ni%d" % len(pos), positions=pos)


def ni_particle(length, cutoff, cut_bottom_layers=0):
    """fcc truncated octahedron (Ni) sitting on its (100) facet; optional flattening."""
    p = truncated_octahedron(length, cutoff)
    z = p.positions[:, 2]
    layers = np.unique(np.round(z, 2))
    if cut_bottom_layers:
        p = p[z > layers[cut_bottom_layers - 1] + 0.1]
    p.positions -= p.positions.min(axis=0)
    p.cell = None
    p.pbc = False
    return p


def classify_sites(p, support_contact=True):
    """Metal coordination and site type for every atom of a particle.

    Types: interface (bottom layer touching the oxide), corner (CN<=6), edge (CN=7),
    face100 (CN=8), face111 (CN=9), bulk (CN>=10, inc. sub-surface).
    """
    nn = A_NI / np.sqrt(2)
    pairs = cKDTree(p.positions).query_pairs(1.2 * nn, output_type="ndarray")
    cn = np.bincount(pairs.ravel(), minlength=len(p))
    z = p.positions[:, 2]
    bottom = z < z.min() + 0.3 if support_contact else np.zeros(len(p), bool)
    t = np.full(len(p), "bulk", dtype=object)
    t[cn == 9] = "face111"
    t[cn == 8] = "face100"
    t[cn == 7] = "edge"
    t[cn <= 6] = "corner"
    t[bottom] = "interface"
    return cn, t


def decorate_cu(p, types, fill=("corner", "edge"), fraction=1.0, seed=2):
    p = p.copy()
    rng = np.random.default_rng(seed)
    idx = np.where(np.isin(types, fill))[0]
    if fraction < 1.0:
        idx = rng.choice(idx, int(round(fraction * len(idx))), replace=False)
    p.symbols[idx] = "Cu"
    return p


def place_on(slab, particle, gap=2.0):
    p = particle.copy()
    top = slab.positions[:, 2].max()
    c = slab.cell[0, :2] / 2 + slab.cell[1, :2] / 2
    p.positions[:, :2] += c - p.positions[:, :2].mean(axis=0)
    p.positions[:, 2] += top + gap - p.positions[:, 2].min()
    sys_ = slab + p
    sys_.cell[2, 2] = max(slab.cell[2, 2], p.positions[:, 2].max() + 8)
    return sys_
