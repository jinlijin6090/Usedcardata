#!/usr/bin/env python
"""Create and monitor all simulation cases.

    python campaign/campaign.py make     [task ...]   create runs/<task>/<case>/config.json
    python campaign/campaign.py list     [task ...]   print case directories (one per line, for job arrays)
    python campaign/campaign.py status   [task ...]   progress table
    python campaign/campaign.py estimate [--sec-per-ps 20]   simulated time and wall-time estimate
    python campaign/campaign.py pending  [task ...]   unfinished cases: "<remaining_ps> <ready 0/1> <dir>"

tasks: task1_flux, task2_anneal, task3_codep   (default: all)
Existing config.json files are never overwritten unless --force is given (protects running cases).
"""
import argparse
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RUNS = os.path.join(ROOT, "runs")
TASKS = ["task1_flux", "task2_anneal", "task3_codep"]


def load():
    return json.load(open(os.path.join(ROOT, "campaign", "campaign.json")))


def fmt(x):
    return f"{x:g}".replace(".", "p")


def cases(camp):
    """Yield (task, case_name, config dict)."""
    base = dict(ffield=camp["ffield"], T=camp["T"], checkpoint_every_ps=camp["checkpoint_every_ps"])
    eq = camp.get("equil_ps", 10)
    t1 = camp["task1_flux"]
    for iv in t1["c_interval_ps"]:
        for s in t1["seeds"]:
            yield "task1_flux", f"i{fmt(iv)}_s{s}", dict(
                base, init_data=camp["cluster"], seed=s,
                equil_ps=eq, deposit={"C": t1["nC"]}, c_interval_ps=iv, post_ps=t1["post_ps"],
                dump_every_ps=max(0.5, iv / 10))
    t2 = camp["task2_anneal"]
    iv = t2["from_c_interval_ps"]
    for T in t2["T"]:
        for s in t2["seeds"]:
            yield "task2_anneal", f"from_i{fmt(iv)}_s{s}_T{fmt(T)}", dict(
                base, T=T, init_data=os.path.join("..", "..", "task1_flux", f"i{fmt(iv)}_s{s}", "final.data"),
                seed=100 + s, equil_ps=0, deposit={}, post_ps=t2["anneal_ps"],
                dump_every_ps=t2["dump_every_ps"], recenter_every_ps=10)
    t3 = camp["task3_codep"]
    for r in t3["O_per_C"]:
        for s in t3["seeds"]:
            nO = int(round(r * t3["nC"]))
            yield "task3_codep", f"OC{fmt(r)}_s{s}", dict(
                base, init_data=camp["cluster"], seed=s,
                equil_ps=eq, deposit={"C": t3["nC"], "O": nO}, c_interval_ps=t3["c_interval_ps"],
                post_ps=t3["post_ps"], dump_every_ps=max(0.5, t3["c_interval_ps"] / 10),
                O_per_C=r)


def total_ps(cfg):
    nC = cfg["deposit"].get("C", 0)
    n = sum(cfg["deposit"].values())
    dep = (nC if nC else n) * cfg.get("c_interval_ps", 0) if n else 0
    return cfg["equil_ps"] + dep + cfg["post_ps"]


def selected(argv_tasks):
    return argv_tasks or TASKS


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=["make", "list", "status", "estimate", "pending"])
    ap.add_argument("tasks", nargs="*", default=[])
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--sec-per-ps", type=float, default=20.0,
                    help="measured wall seconds per ps (see progress.log); ~11-25 s/ps on 4 threads")
    a = ap.parse_args()
    bad = [t for t in a.tasks if t not in TASKS]
    if bad:
        sys.exit(f"unknown task(s) {bad}; choose from {TASKS}")
    camp = load()
    tasks = selected(a.tasks)
    allc = [c for c in cases(camp) if c[0] in tasks]

    if a.cmd == "make":
        n_new = 0
        for task, name, cfg in allc:
            d = os.path.join(RUNS, task, name)
            os.makedirs(d, exist_ok=True)
            f = os.path.join(d, "config.json")
            if os.path.exists(f) and not a.force:
                continue
            cfg = dict(cfg, name=name, task=task)
            json.dump(cfg, open(f, "w"), indent=1)
            n_new += 1
        print(f"{n_new} new config(s) written under {RUNS}  ({len(allc)} cases in selected tasks)")

    elif a.cmd == "list":
        for task, name, _ in allc:
            print(os.path.join(RUNS, task, name))

    elif a.cmd == "status":
        print(f"{'case':38s} {'phase':8s} {'t (ps)':>9s} {'of':>8s} {'%':>5s}  last progress line")
        for task, name, cfg in allc:
            d = os.path.join(RUNS, task, name)
            st = os.path.join(d, "state.json")
            tot = total_ps(cfg)
            if os.path.exists(st):
                s = json.load(open(st))
                t = s["step"] * cfg.get("dt_fs", 0.25) / 1000
                phase = s["phase"]
            else:
                t, phase = 0.0, ("not made" if not os.path.exists(os.path.join(d, "config.json"))
                                 else "running?" if os.path.exists(os.path.join(d, "progress.log"))
                                 else "waiting")
            last = ""
            pl = os.path.join(d, "progress.log")
            if os.path.exists(pl):
                lines = open(pl).read().strip().splitlines()
                last = lines[-1][:90] if lines else ""
            print(f"{task+'/'+name:38s} {phase:8s} {t:9.1f} {tot:8.0f} {100*t/tot:5.1f}  {last}")

    elif a.cmd == "pending":
        for task, name, cfg in allc:
            d = os.path.join(RUNS, task, name)
            if not os.path.exists(os.path.join(d, "config.json")):
                continue
            cfg = json.load(open(os.path.join(d, "config.json")))
            st = os.path.join(d, "state.json")
            t = 0.0
            if os.path.exists(st):
                s = json.load(open(st))
                if s.get("phase") == "done":
                    continue
                t = s["step"] * cfg.get("dt_fs", 0.25) / 1000
            init = cfg["init_data"]
            ready = int(os.path.exists(os.path.join(d, init)) or os.path.exists(os.path.join(ROOT, init))
                        or os.path.exists(os.path.join(d, "ckpt.data")))
            print(f"{total_ps(cfg) - t:.1f} {ready} {d}")

    elif a.cmd == "estimate":
        print(f"assuming {a.sec_per_ps} wall s per simulated ps per case (4 OpenMP threads)\n")
        print(f"{'task':14s} {'cases':>5s} {'ps/case (max)':>14s} {'total ps':>10s} {'wall h (longest case)':>22s} {'core-h (4 thr)':>15s}")
        for task in tasks:
            cs = [total_ps(c[2]) for c in allc if c[0] == task]
            if not cs:
                continue
            tot = sum(cs)
            print(f"{task:14s} {len(cs):5d} {max(cs):14.0f} {tot:10.0f} "
                  f"{max(cs)*a.sec_per_ps/3600:22.1f} {tot*a.sec_per_ps/3600*4:15.0f}")


if __name__ == "__main__":
    main()
