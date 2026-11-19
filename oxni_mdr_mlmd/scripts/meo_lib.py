"""Model builders for the OxNi dry-reforming catalyst (He et al., ACS Catal. 2021, 11, 12409).

Paper picture: ~15 nm spherical-like Ni particles, covered by a 1-2 nm multielement oxide
(MEO: Ni2+, Mg, Al, Si, O; rock-salt-like, d = 2.2 A), sitting on a MgO-Al2O3 matrix
from a calcined Ni-Mg-Al hydrotalcite. The MEO is porous, so some surface Ni stays exposed.

Model hierarchy built here (see PLAN.md):
  slab   : Ni(211)/(111)/(100) slab + 1-1.5 nm MEO film with a pore that exposes Ni,
           + CH4/CO2 gas or pre-adsorbed C*/O*/H*. Bare Ni slab as the reference catalyst.
  particle: fcc Ni sphere (flat (100) contact facet) on a Mg(Al)O rock-salt support,
           optionally wrapped in a porous MEO shell. Bare particle as the reference.

The MEO starts as a random rock-salt solid solution with cation/anion vacancies chosen so the
oxide is formally charge neutral (Ni2+, Mg2+, Al3+, Si4+, O2-). This is only an initial guess:
anneal it with the ML potential (protocols/meo_anneal.json) before production runs.

Tag convention (atoms.get_tags()):
  0 Ni, 1 Ni frozen, 2 MEO, 3 support, 4 support frozen, 5 gas molecule, 6 pre-adsorbed atom
"""
import numpy as np
from ase import Atoms
from ase.build import fcc100, fcc111, fcc211
from scipy.spatial.transform import Rotation

A_NI = 3.524
A_RS = 4.20  # MEO rock-salt lattice parameter (MgO 4.21, NiO 4.18); paper: d = 2.2 A fringes
A_MGO = 4.212

TAG = dict(ni=0, ni_fixed=1, meo=2, support=3, support_fixed=4, gas=5, ads=6)
CHARGE = {"Mg": 2, "Ni": 2, "Al": 3, "Si": 4, "O": -2}

# MEO cation composition (fractions). The paper only says the overlayer holds Ni, Mg, Al and Si
# in a high-entropy-like rock-salt phase; equimolar is the working hypothesis. Calibrate with
# EDS/XPS ratios and treat it as a sensitivity parameter.
MEO_EQUIMOLAR = {"Mg": 0.25, "Ni": 0.25, "Al": 0.25, "Si": 0.25}
# Support: Mg(Al)O periclase-type mixed oxide from the hydrotalcite (Mg:Al = 12:5 in the recipe).
SUPPORT_MGALO = {"Mg": 12 / 17, "Al": 5 / 17}

COLORS = {  # roughly the palette of Fig. 1i of the paper
    "Ni": (0.70, 0.13, 0.15), "Mg": (0.95, 0.76, 0.28), "Al": (0.45, 0.75, 0.30),
    "Si": (0.35, 0.60, 0.90), "O": (0.88, 0.88, 0.88), "C": (0.20, 0.20, 0.20),
    "H": (1.00, 1.00, 1.00),
}


# ----------------------------------------------------------------------------- Ni slabs
def ni_slab(facet="211", nx=9, ny=12, layers=4, fix_bottom=3.0, vacuum=12.0, a=A_NI):
    """Orthogonal Ni slab; atoms within `fix_bottom` A of the bottom are tagged frozen.

    facet '211' (step + (111) terrace), '111' (ny must be even) or '100'.
    For '211' (ASE fcc211) nx must be a multiple of 3 and `layers` counts 3-row step units
    (~2.2 A each), so layers=4 is ~7.9 A thick; for '111'/'100' it counts atomic layers.
    """
    if facet == "111":
        s = fcc111("Ni", (nx, ny, layers), a=a, orthogonal=True)
    elif facet == "100":
        s = fcc100("Ni", (nx, ny, layers), a=a)
    elif facet == "211":
        s = fcc211("Ni", (nx, ny, layers), a=a, orthogonal=True)
    else:
        raise ValueError(facet)
    s.positions[:, 2] -= s.positions[:, 2].min()
    s.set_tags(np.where(s.positions[:, 2] < fix_bottom, TAG["ni_fixed"], TAG["ni"]))
    s.cell[2, 2] = s.positions[:, 2].max() + vacuum
    s.pbc = (True, True, False)
    return s


def _unit_lengths(facet, a=A_NI):
    base = {"111": (1, 2), "100": (1, 1), "211": (3, 1)}[facet]
    s = ni_slab(facet, base[0], base[1], layers=3, fix_bottom=0.0, a=a)
    return s.cell[0, 0], s.cell[1, 1], base


def match_film_to_slab(facet, lmin=18.0, lmax=32.0, a_rs=A_RS, a=A_NI):
    """Search slab repeats (nx, ny) and film repeats (px, py) with the smallest film strain.

    The rock-salt (100) film uses its square surface cell (a_rs / sqrt 2). Returns a list of
    dicts sorted by max |strain| (strain = slab length / film length - 1).
    """
    ux, uy, (bx, by) = _unit_lengths(facet, a)
    a_s = a_rs / np.sqrt(2)
    out = []
    for kx in range(1, int(lmax / ux) + 1):
        for ky in range(1, int(lmax / uy) + 1):
            lx, ly = kx * ux, ky * uy
            if not (lmin <= lx <= lmax and lmin <= ly <= lmax):
                continue
            px, py = max(1, round(lx / a_s)), max(1, round(ly / a_s))
            sx, sy = lx / (px * a_s) - 1, ly / (py * a_s) - 1
            out.append(dict(nx=kx * bx, ny=ky * by, px=px, py=py, lx=lx, ly=ly,
                            strain_x=sx, strain_y=sy, max_strain=max(abs(sx), abs(sy))))
    return sorted(out, key=lambda d: (round(d["max_strain"], 3), d["lx"] * d["ly"]))


# ------------------------------------------------------------------ oxide decoration
def _species_counts(composition, n):
    names = list(composition)
    f = np.array([composition[k] for k in names], float)
    f /= f.sum()
    raw = f * n
    cnt = np.floor(raw).astype(int)
    for i in np.argsort(-(raw - cnt))[: n - cnt.sum()]:
        cnt[i] += 1
    return dict(zip(names, cnt))


def _decorate(cat_pos, an_pos, composition, o_deficit, rng, tag):
    """Fill rock-salt cation/anion sites so that the oxide is formally neutral.

    O vacancies: fraction `o_deficit` of anion sites (reduced oxide, e.g. after H2 at 1073 K).
    Cation vacancies follow from charge balance (Al3+, Si4+ on a rock-salt lattice).
    """
    n_cs, n_as = len(cat_pos), len(an_pos)
    n_o = n_as - int(round(o_deficit * n_as))
    q = {k: CHARGE[k] for k in composition}
    f = np.array(list(composition.values()), float)
    q_avg = np.dot(f / f.sum(), list(q.values()))
    n_c = min(n_cs, int(round(2 * n_o / q_avg)))
    cnt = _species_counts(composition, n_c)

    def total():
        return sum(cnt[k] * q[k] for k in cnt) - 2 * n_o

    # make the total charge exactly zero: odd -> swap a 2+/3+ pair, even -> adjust O count
    for _ in range(100):
        Q = total()
        if Q == 0:
            break
        if Q % 2:
            if "Al" in cnt and "Mg" in cnt and Q > 0 and cnt["Al"] > 0:
                cnt["Al"] -= 1; cnt["Mg"] += 1
            elif "Al" in cnt and "Mg" in cnt and Q < 0 and cnt["Mg"] > 0:
                cnt["Mg"] -= 1; cnt["Al"] += 1
            else:
                raise ValueError("cannot balance an odd charge without both Mg and Al")
        elif Q > 0:
            if n_o + Q // 2 <= n_as:
                n_o += Q // 2
            else:
                k = "Mg" if cnt.get("Mg", 0) else max(cnt, key=cnt.get)
                cnt[k] -= 1
        else:
            n_o += Q // 2  # Q < 0: remove O
    n_c = sum(cnt.values())
    if n_c > n_cs or n_o > n_as:
        raise ValueError("composition does not fit on the lattice")
    ci = rng.choice(n_cs, n_c, replace=False)
    ai = rng.choice(n_as, n_o, replace=False)
    sym = np.concatenate([np.repeat(list(cnt), list(cnt.values()))])
    rng.shuffle(sym)
    at = Atoms(list(sym) + ["O"] * n_o, positions=np.vstack([cat_pos[ci], an_pos[ai]]))
    at.set_tags(np.full(len(at), tag))
    return at


def rocksalt_film(px, py, n_layers, composition=MEO_EQUIMOLAR, a=A_RS, o_deficit=0.0,
                  seed=0, tag=TAG["meo"]):
    """Random rock-salt (100) film, px x py square surface cells (a/sqrt2), n_layers thick."""
    rng = np.random.default_rng(seed)
    a_s = a / np.sqrt(2)
    i, j = np.meshgrid(np.arange(px), np.arange(py), indexing="ij")
    i, j = i.ravel(), j.ravel()
    cat, an = [], []
    for k in range(n_layers):
        off = 0.5 * (k % 2)
        z = np.full(i.size, k * a / 2)
        cat.append(np.c_[(i + off) * a_s, (j + off) * a_s, z])
        an.append(np.c_[(i + 0.5 - off) * a_s, (j + 0.5 - off) * a_s, z])
    f = _decorate(np.vstack(cat), np.vstack(an), composition, o_deficit, rng, tag)
    f.cell = [px * a_s, py * a_s, (n_layers - 1) * a / 2 + 1.0]
    f.pbc = (True, True, False)
    return f


def region_charge(atoms, tag):
    m = atoms.get_tags() == tag
    return int(sum(CHARGE.get(s, 0) for s in atoms.symbols[m]))


def _xy_dist(atoms, center_xy):
    L = np.diag(atoms.cell)[:2]
    d = atoms.positions[:, :2] - np.asarray(center_xy)
    d -= np.round(d / L) * L
    return np.linalg.norm(d, axis=1)


def neutralize(atoms, tag, order=None):
    """Remove O / divalent cations (or swap Mg<->Al) of region `tag` until it is neutral.

    `order(atoms, idx)` returns idx sorted by removal priority (default: as given).
    """
    s = atoms.copy()
    for _ in range(10000):
        Q = region_charge(s, tag)
        if Q == 0:
            return s
        idx = np.where(s.get_tags() == tag)[0]
        if order is not None:
            idx = order(s, idx)
        sym = s.symbols
        if Q <= -2:
            del s[[i for i in idx if sym[i] == "O"][0]]
        elif Q >= 2:
            c = [i for i in idx if sym[i] in ("Mg", "Ni")]
            if c:
                del s[c[0]]
            else:
                s.symbols[[i for i in idx if sym[i] in ("Al", "Si")][0]] = "Mg" if Q == 2 else "Al"
        elif Q == -1:
            c = [i for i in idx if sym[i] == "Mg"] or [i for i in idx if sym[i] == "Al"]
            s.symbols[c[0]] = "Al" if sym[c[0]] == "Mg" else "Si"
        else:  # Q == 1
            c = [i for i in idx if sym[i] == "Al"] or [i for i in idx if sym[i] == "Si"]
            s.symbols[c[0]] = "Mg" if sym[c[0]] == "Al" else "Al"
    raise RuntimeError("neutralize did not converge")


def stack_film(substrate, film, gap=2.1, vacuum=12.0):
    """Strain the film onto the substrate cell and put it `gap` A above the substrate top."""
    f = film.copy()
    sx = substrate.cell[0, 0] / f.cell[0, 0]
    sy = substrate.cell[1, 1] / f.cell[1, 1]
    f.positions[:, 0] *= sx
    f.positions[:, 1] *= sy
    f.positions[:, 2] += substrate.positions[:, 2].max() + gap - f.positions[:, 2].min()
    out = substrate.copy()
    out += f
    out.cell[2, 2] = out.positions[:, 2].max() + vacuum
    out.pbc = (True, True, False)
    return out, (sx - 1, sy - 1)


def carve_pore(system, radius, center_xy=None, tag=TAG["meo"]):
    """Cylindrical hole through the MEO film (exposes Ni). Re-neutralises from the rim outwards."""
    s = system.copy()
    c = np.diag(s.cell)[:2] / 2 if center_xy is None else np.asarray(center_xy)
    rm = (s.get_tags() == tag) & (_xy_dist(s, c) < radius)
    del s[np.where(rm)[0]]
    s = neutralize(s, tag, order=lambda a, idx: idx[np.argsort(_xy_dist(a, c)[idx])])
    s.info["pore_center_xy"] = tuple(map(float, c))
    s.info["pore_radius"] = float(radius)
    return s


# ------------------------------------------------------------------- gas & adsorbates
def ch4():
    v = np.array([[1, 1, 1], [-1, -1, 1], [-1, 1, -1], [1, -1, -1]]) / np.sqrt(3)
    return Atoms("CH4", positions=np.vstack([[0, 0, 0], 1.09 * v]))


def co2():
    return Atoms("CO2", positions=[[0, 0, 0], [0, 0, 1.16], [0, 0, -1.16]])


def _min_dist_xy_periodic(atoms, pos):
    L = np.diag(atoms.cell)[:2]
    if len(atoms) == 0:
        return np.inf
    d = atoms.positions[None, :, :] - pos[:, None, :]
    d[..., :2] -= np.round(d[..., :2] / L) * L
    return np.linalg.norm(d, axis=-1).min()


def add_gas(system, n_ch4, n_co2, z_gap=3.0, height=20.0, min_dist=3.0, seed=0, wall_pad=1.5):
    """Random CH4/CO2 molecules in a layer `height` A thick, starting `z_gap` above the top.

    Sets atoms.info['wall_z'] (upper reflecting wall used by run_md.py). Note: at 1 atm and
    873 K a 25x25x20 A box holds ~0.1 molecule; any practical loading is ~10^2-10^3 bar,
    so read rates per collision, not per unit time.
    """
    rng = np.random.default_rng(seed)
    s = system.copy()
    zlo = s.positions[:, 2].max() + z_gap
    zhi = zlo + height
    L = np.diag(s.cell)[:2]
    kinds = ["CH4"] * n_ch4 + ["CO2"] * n_co2
    rng.shuffle(kinds)
    for k in kinds:
        for _ in range(5000):
            m = ch4() if k == "CH4" else co2()
            m.positions = Rotation.random(random_state=rng.integers(1 << 31)).apply(m.positions)
            com = np.r_[rng.random(2) * L, zlo + 1.2 + rng.random() * (height - 2.4)]
            m.positions += com
            if _min_dist_xy_periodic(s, m.positions) > min_dist:
                m.set_tags(np.full(len(m), TAG["gas"]))
                s += m
                break
        else:
            raise RuntimeError("could not place gas molecule; enlarge `height`")
    s.info["wall_z"] = float(zhi + wall_pad)
    s.cell[2, 2] = max(s.cell[2, 2], zhi + wall_pad + 4.0)
    return s


def add_adatoms(system, counts, center_xy=None, radius=None, bond=None, min_dist=2.4,
                min_ni=3, seed=0):
    """Put C*/O*/H* on exposed Ni (inside the disc center_xy/radius if given).

    Each adatom is lowered onto the Ni surface until its nearest Ni is `bond` A away and
    must have >= `min_ni` Ni within bond + 0.3 A (hollow/bridge-like sites).
    """
    bond = bond or {"C": 1.80, "O": 1.88, "H": 1.70}
    rng = np.random.default_rng(seed)
    s = system.copy()
    L = np.diag(s.cell)[:2]
    c = np.diag(s.cell)[:2] / 2 if center_xy is None else np.asarray(center_xy)
    ni_pos = s.positions[np.isin(s.get_tags(), [TAG["ni"], TAG["ni_fixed"]])]
    for el, n in counts.items():
        b = bond[el]
        for _ in range(n):
            for _try in range(20000):
                if radius is None:
                    xy = rng.random(2) * L
                else:
                    r, t = radius * np.sqrt(rng.random()), 2 * np.pi * rng.random()
                    xy = (c + r * np.array([np.cos(t), np.sin(t)])) % L
                d = ni_pos[:, :2] - xy
                d -= np.round(d / L) * L
                dxy = np.linalg.norm(d, axis=1)
                under = dxy < b
                if not under.any():
                    continue
                z = (ni_pos[under, 2] + np.sqrt(b ** 2 - dxy[under] ** 2)).max()
                pos = np.r_[xy, z]
                dist_ni = np.sqrt(dxy ** 2 + (ni_pos[:, 2] - z) ** 2)
                if (dist_ni < b + 0.3).sum() < min_ni:
                    continue
                other = s[~np.isin(s.get_tags(), [TAG["ni"], TAG["ni_fixed"]])]
                other.cell = s.cell
                if _min_dist_xy_periodic(other, pos[None]) > min_dist:
                    s += Atoms(el, positions=[pos], tags=[TAG["ads"]])
                    break
            else:
                raise RuntimeError("could not place adatom %s" % el)
    return s


# ------------------------------------------------------------------- particle models
def _cubic_grid(rmax, a):
    n = int(np.ceil(rmax / (a / 2))) + 1
    r = np.arange(-n, n + 1)
    i, j, k = [x.ravel() for x in np.meshgrid(r, r, r, indexing="ij")]
    return np.c_[i, j, k]


def ni_sphere(diameter, cut=0.15, a=A_NI):
    """Spherical fcc Ni (high-index, stepped surface as in the paper) with a flat (100)
    contact facet: the bottom `cut` fraction of the height is removed. Centered at 0."""
    R = diameter / 2
    g = _cubic_grid(R, a)
    g = g[(g.sum(axis=1) % 2) == 0]
    p = g * (a / 2)
    p = p[(np.linalg.norm(p, axis=1) <= R) & (p[:, 2] >= -R + cut * diameter)]
    at = Atoms("Ni%d" % len(p), positions=p, tags=[TAG["ni"]] * len(p))
    at.info["R"] = R
    return at


def _pore_axes(n, rng, zmin=0.15, min_angle=50.0):
    axes = []
    for _ in range(20000):
        if len(axes) == n:
            break
        v = rng.normal(size=3)
        v /= np.linalg.norm(v)
        if v[2] < zmin:
            continue
        if all(np.degrees(np.arccos(np.clip(v @ w, -1, 1))) > min_angle for w in axes):
            axes.append(v)
    return np.array(axes)


def meo_shell(r_in, thickness, composition=MEO_EQUIMOLAR, a=A_RS, z_floor=-np.inf,
              n_pores=6, pore_radius=4.5, o_deficit=0.0, seed=0):
    """Porous rock-salt MEO shell (cube-on-cube with Ni) between r_in and r_in+thickness,
    above z_floor (relative to the sphere center). Pores are cylinders along random
    upper-hemisphere directions."""
    rng = np.random.default_rng(seed)
    g = _cubic_grid(r_in + thickness, a)
    p = g * (a / 2)
    r = np.linalg.norm(p, axis=1)
    keep = (r >= r_in) & (r <= r_in + thickness) & (p[:, 2] >= z_floor)
    axes = _pore_axes(n_pores, rng) if n_pores else np.zeros((0, 3))
    for v in axes:
        t = p @ v
        perp = np.linalg.norm(p - np.outer(t, v), axis=1)
        keep &= ~((t > 0) & (perp < pore_radius))
    even = (g.sum(axis=1) % 2) == 0
    sh = _decorate(p[keep & even], p[keep & ~even], composition, o_deficit, rng, TAG["meo"])
    sh.info["pore_axes"] = axes.tolist()
    return sh


def rocksalt_support(px, py, n_layers=4, n_fixed=2, composition=SUPPORT_MGALO, a=A_MGO,
                     seed=0):
    s = rocksalt_film(px, py, n_layers, composition, a=a, seed=seed, tag=TAG["support"])
    z = s.positions[:, 2]
    s.set_tags(np.where(z < (n_fixed - 0.5) * a / 2, TAG["support_fixed"], TAG["support"]))
    return s


def particles_on_support(particles, support_cells, support_layers=4, gap=2.1, shell_gap=2.0,
                         meo_composition=MEO_EQUIMOLAR, o_deficit=0.0, vacuum=15.0, seed=0):
    """Ni spheres (optionally MEO-wrapped) on a Mg(Al)O(100) support.

    particles: list of dicts {xy: (fx, fy) fractional, d: Ni diameter (A),
               shell: MEO thickness (A) or 0, n_pores, pore_radius}
    """
    sup = rocksalt_support(*support_cells, n_layers=support_layers, seed=seed)
    top = sup.positions[:, 2].max()
    L = np.diag(sup.cell)[:2]
    parts, shells = [], []
    for k, spec in enumerate(particles):
        ni = ni_sphere(spec["d"])
        c = np.r_[np.asarray(spec["xy"]) * L, top + gap - ni.positions[:, 2].min()]
        ni.positions += c
        parts.append(ni)
        if spec.get("shell"):
            sh = meo_shell(ni.info["R"] + shell_gap, spec["shell"], meo_composition,
                           z_floor=top + 1.6 - c[2], n_pores=spec.get("n_pores", 6),
                           pore_radius=spec.get("pore_radius", 4.5), o_deficit=o_deficit,
                           seed=seed + 11 * k + 1)
            sh.positions += c
            shells.append(sh)
    system = sup.copy()
    for p in parts:
        system += p
    if shells:
        meo = shells[0]
        for sh in shells[1:]:
            meo += sh
        # drop shell sites that clash with Ni/support or with another shell
        from scipy.spatial import cKDTree
        tree = cKDTree(system.positions)
        clash = np.array([len(x) > 0 for x in tree.query_ball_point(meo.positions, 1.8)])
        meo = meo[~clash]
        pairs = cKDTree(meo.positions).query_pairs(1.6, output_type="ndarray")
        if len(pairs):
            del meo[np.unique(pairs[:, 1])]
        system += meo
        system = neutralize(system, TAG["meo"],
                            order=lambda a, idx: idx[np.argsort(-a.positions[idx, 2])])
    system.cell[2, 2] = system.positions[:, 2].max() + vacuum
    system.pbc = (True, True, False)
    return system


# ------------------------------------------------------------------------- reporting
def summary(atoms):
    names = {v: k for k, v in TAG.items()}
    tags = atoms.get_tags()
    lines = ["%d atoms, cell %s" % (len(atoms), np.round(atoms.cell.lengths(), 2))]
    for t in np.unique(tags):
        m = tags == t
        f = atoms[m].get_chemical_formula(mode="metal")
        q = region_charge(atoms, t) if t == TAG["meo"] else None
        if t == TAG["support"]:
            q = region_charge(atoms, t) + region_charge(atoms, TAG["support_fixed"])
        lines.append("  tag %d %-13s %5d  %s%s" % (t, names[t], m.sum(), f,
                                                  "" if q is None else "  (formal q=%+d)" % q))
    return "\n".join(lines)


def colors_for(atoms):
    return [COLORS.get(s, (0.6, 0.6, 0.6)) for s in atoms.symbols]
