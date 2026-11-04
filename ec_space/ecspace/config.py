"""Configuration: INI file + command-line overrides (--set section.key=value)."""
from __future__ import annotations

import configparser
import os
from dataclasses import dataclass, field
from pathlib import Path


def _list(s: str) -> list[str]:
    return [x.strip() for x in s.replace("\n", ",").split(",") if x.strip()]


def _bool(s: str) -> bool:
    return s.strip().lower() in ("1", "yes", "true", "on", "y")


def detect_workers() -> int:
    for var in ("SLURM_CPUS_PER_TASK", "SLURM_CPUS_ON_NODE", "PBS_NUM_PPN", "NCPUS", "LSB_DJOB_NUMPROC", "NSLOTS", "OMP_NUM_THREADS_EC"):
        v = os.environ.get(var)
        if v and v.split("(")[0].isdigit() and int(v.split("(")[0]) > 0:
            return int(v.split("(")[0])
    try:
        return len(os.sched_getaffinity(0))
    except AttributeError:
        return os.cpu_count() or 1


@dataclass
class Config:
    max_heavy: int = 12
    max_C: int = 12
    max_O: int = 12
    max_H: int | None = None
    allow_peroxide: bool = True
    strain_filter: str = "valence"
    special: list[str] = field(default_factory=lambda: ["[C-]#[O+]"])
    start: str = "O=C1OCCO1"
    radicals: bool = True
    anions: bool = True
    cation_centers: list[str] = field(default_factory=lambda: ["C"])
    protonation: bool = True
    radical_ions: bool = True
    pairs: list[str] = field(default_factory=lambda: ["distonic", "dianion"])
    carbenes: bool = True
    li_max: int = 2
    charges: list[int] = field(default_factory=lambda: [-2, -1, 0, 1])
    families: list[str] = field(default_factory=list)
    derive_max_heavy: int = 9
    count_sample: float = 0.002
    gzip_level: int = 3
    out: Path = Path("ec_space_out")
    workers: int = 1
    shards: str = "auto"
    dry_run: bool = False
    viewer_max_heavy: int = 7
    max_per_cell: int = 200000
    max_species: int = 2000000
    source: str = ""


def load(path: str | None, overrides: list[str] | None = None) -> Config:
    cp = configparser.ConfigParser(inline_comment_prefixes=(";",), interpolation=None)
    base = Path.cwd()
    if path:
        p = Path(path).resolve()
        cp.read(p, encoding="utf-8")
        base = p.parent
    for ov in overrides or []:
        k, v = ov.split("=", 1)
        sec, key = k.split(".", 1)
        if not cp.has_section(sec):
            cp.add_section(sec)
        cp.set(sec, key, v)
    g = lambda s, k, d="": cp.get(s, k, fallback=d).strip()  # noqa: E731
    c = Config()
    c.max_heavy = int(g("space", "max_heavy", "12"))
    c.max_C = int(g("space", "max_C", str(c.max_heavy)) or c.max_heavy)
    c.max_O = int(g("space", "max_O", str(c.max_heavy)) or c.max_heavy)
    h = g("space", "max_H", "")
    c.max_H = int(h) if h else None
    c.allow_peroxide = _bool(g("space", "allow_peroxide", "yes"))
    c.strain_filter = g("space", "strain_filter", "valence") or "valence"
    c.special = _list(g("space", "special", "[C-]#[O+]"))
    c.start = g("space", "start", "O=C1OCCO1") or "O=C1OCCO1"
    c.radicals = _bool(g("species", "radicals", "yes"))
    c.anions = _bool(g("species", "anions", "yes"))
    c.cation_centers = [x.upper() for x in _list(g("species", "cation_centers", "C"))]
    c.protonation = _bool(g("species", "protonation", "yes"))
    c.radical_ions = _bool(g("species", "radical_ions", "yes"))
    c.pairs = _list(g("species", "pairs", "distonic, dianion"))
    c.carbenes = _bool(g("species", "carbenes", "yes"))
    c.li_max = int(g("species", "li_max", "2") or 0)
    c.charges = [int(x) for x in _list(g("species", "charges", "-2,-1,0,1"))]
    fam = g("reactions", "families", "")
    c.families = _list(fam) if fam else ["h_loss", "deprotonation", "hydride_loss", "protonation", "electron_capture",
                                         "electron_loss", "site_electron_transfer", "li_binding", "homolysis", "heterolysis",
                                         "beta_scission", "radical_anion_fragmentation", "radical_cation_fragmentation", "pair_steps"]
    c.derive_max_heavy = int(g("write", "derive_max_heavy", "9"))
    c.count_sample = float(g("write", "count_sample", "0.002") or 0)
    c.gzip_level = int(g("write", "gzip_level", "3"))
    out = Path(g("run", "out", "./ec_space_out") or "./ec_space_out")
    c.out = out if out.is_absolute() else (base / out).resolve()
    w = g("run", "workers", "auto")
    c.workers = detect_workers() if w in ("", "auto") else int(w)
    c.shards = g("run", "shards", "auto") or "auto"
    c.dry_run = _bool(g("run", "dry_run", "no"))
    c.viewer_max_heavy = int(g("viewer", "max_heavy", "7"))
    c.max_per_cell = int(g("viewer", "max_per_cell", "200000"))
    c.max_species = int(g("viewer", "max_species", "2000000"))
    c.source = str(path or "")
    return c
