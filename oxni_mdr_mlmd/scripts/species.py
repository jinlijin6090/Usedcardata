"""Graph-based species / coke / structure analysis for C-H-O chemistry on Ni + oxide.

Bonds between C, H and O use element-pair cutoffs. Lattice O (tags MEO/support) is only
counted when it is bonded to C or H (carbonate, OH, CO made from lattice O = Mars-van
Krevelen marker). A fragment is "adsorbed" (*) when any of its atoms is within the
surface cutoff of a cation (Ni, Mg, Al, Si) or when it contains lattice O.
"""
from collections import Counter

import numpy as np
from ase.neighborlist import neighbor_list
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

BOND = {("C", "H"): 1.35, ("C", "O"): 1.60, ("H", "O"): 1.25, ("H", "H"): 1.00,
        ("C", "C"): 1.85, ("O", "O"): 1.55}
SURF = {"C": 2.30, "O": 2.40, "H": 2.00}
CATIONS = ("Ni", "Mg", "Al", "Si")
LATTICE_TAGS = (2, 3, 4)
NI_TAGS = (0, 1)
RMAX = 2.45

KEY_SPECIES = ["CH4(g)", "CO2(g)", "CO(g)", "H2(g)", "H2O(g)",
               "CH4*", "CH3*", "CH2*", "CH*", "C*", "H*", "O*", "OH*", "CO*", "CO2*",
               "CHO*", "COH*", "CHO2*", "CO3*", "H2O*"]


def _pair_cut(a, b):
    return BOND.get((a, b)) or BOND.get((b, a))


def formula(sym):
    c = Counter(sym)
    out = ""
    for el in ("C", "H", "O"):
        n = c.get(el, 0)
        out += "" if n == 0 else (el if n == 1 else "%s%d" % (el, n))
    return {"HO": "OH", "CHO2": "CHO2", "CH2O": "CH2O"}.get(out, out)


def analyse(atoms, nl=None):
    """Return (species Counter, metrics dict) for one frame."""
    sym = np.array(atoms.get_chemical_symbols())
    tags = atoms.get_tags()
    i, j, d = nl if nl is not None else neighbor_list("ijd", atoms, RMAX)
    lat_o = (sym == "O") & np.isin(tags, LATTICE_TAGS)
    light = np.isin(sym, ("C", "H", "O"))
    is_cat = np.isin(sym, CATIONS)

    # light-light bonds (no lattice-O / lattice-O pairs)
    m = light[i] & light[j] & ~(lat_o[i] & lat_o[j])
    li, lj, ld = i[m], j[m], d[m]
    cut = np.array([_pair_cut(a, b) for a, b in zip(sym[li], sym[lj])])
    b = ld < cut
    bi, bj = li[b], lj[b]

    # surface contact (adsorbate light atom -> cation)
    m2 = light[i] & is_cat[j]
    si, sd = i[m2], d[m2]
    scut = np.array([SURF[s] for s in sym[si]]) if len(si) else np.zeros(0)
    on_surf = np.zeros(len(atoms), bool)
    on_surf[si[sd < scut]] = True

    nodes = np.where(light & ~lat_o)[0]
    # lattice O enters the graph only when bonded to C or H
    lat_bonded = np.unique(np.r_[bi[lat_o[bi]], bj[lat_o[bj]]])
    nodes = np.union1d(nodes, lat_bonded)
    pos = -np.ones(len(atoms), int)
    pos[nodes] = np.arange(len(nodes))
    keep = (pos[bi] >= 0) & (pos[bj] >= 0)
    g = coo_matrix((np.ones(keep.sum()), (pos[bi[keep]], pos[bj[keep]])),
                   shape=(len(nodes), len(nodes)))
    ncomp, lab = connected_components(g, directed=False)

    species = Counter()
    o_lat_gas = o_lat_c = h_on_oxide = 0
    for c in range(ncomp):
        idx = nodes[lab == c]
        s = sym[idx]
        if len(idx) == 1 and lat_o[idx[0]]:
            continue
        ads = on_surf[idx].any() or lat_o[idx].any()
        name = formula(s) + ("*" if ads else "(g)")
        species[name] += 1
        nlat = int(lat_o[idx].sum())
        if nlat and "C" in s:
            o_lat_c += nlat
            if not ads:
                o_lat_gas += nlat
        if nlat and "C" not in s:
            h_on_oxide += int((s == "H").sum())

    # coke: C-C clusters, dissolved (interstitial) C
    cc = (sym[bi] == "C") & (sym[bj] == "C")
    cidx = np.where(sym == "C")[0]
    largest = n_in_clusters = 0
    if len(cidx):
        cp = -np.ones(len(atoms), int)
        cp[cidx] = np.arange(len(cidx))
        gc = coo_matrix((np.ones(cc.sum()), (cp[bi[cc]], cp[bj[cc]])),
                        shape=(len(cidx), len(cidx)))
        _, cl = connected_components(gc, directed=False)
        sizes = np.bincount(cl)
        largest = int(sizes.max())
        n_in_clusters = int(sizes[sizes >= 2].sum())
    ni = (sym == "Ni")
    m3 = (sym[i] == "C") & ni[j] & (d < 2.15)
    n_ni_c = np.bincount(i[m3], minlength=len(atoms))
    metrics = dict(
        n_C=int(len(cidx)),
        C_in_Cn=n_in_clusters,
        largest_C_cluster=largest,
        C_dissolved=int((n_ni_c[cidx] >= 6).sum()) if len(cidx) else 0,
        O_lattice_in_C_species=o_lat_c,
        O_lattice_in_gas=o_lat_gas,
        H_on_oxide=h_on_oxide,
    )
    return species, metrics


def ni_metrics(atoms, ni_cut=3.0):
    """Ni structure: Rg / mean CN of the metallic Ni (tags 0,1), detached Ni (CN<=2, sintering
    and Ostwald-ripening precursors), metal Ni that went oxidic (>=3 O within 2.3 A, <=3 Ni),
    and MEO Ni2+ that joined the metal (>=6 metal-Ni neighbours, <=1 O; exsolution marker)."""
    sym = np.array(atoms.get_chemical_symbols())
    tags = atoms.get_tags()
    metal = (sym == "Ni") & np.isin(tags, NI_TAGS)
    meo_ni = (sym == "Ni") & (tags == 2)
    if not metal.any():
        return {}
    i, j = neighbor_list("ij", atoms, ni_cut)
    m = metal[j]
    cn_metal = np.bincount(i[m], minlength=len(atoms))
    i2, j2 = neighbor_list("ij", atoms, 2.3)
    mo = (sym[j2] == "O") & np.isin(tags[j2], LATTICE_TAGS + (5, 6))
    n_o = np.bincount(i2[mo], minlength=len(atoms))
    idx = np.where(metal)[0]
    p = atoms.positions[idx]
    rg = float(np.sqrt(((p - p.mean(axis=0)) ** 2).sum(axis=1).mean()))
    return dict(
        Ni_Rg=rg,
        Ni_meanCN=float(cn_metal[idx].mean()),
        Ni_detached=int((cn_metal[idx] <= 2).sum()),
        Ni_oxidic=int(((n_o[idx] >= 3) & (cn_metal[idx] <= 3)).sum()),
        Ni_meo_to_metal=int(((cn_metal >= 6) & (n_o <= 1) & meo_ni).sum()),
    )
