"""Cross-platform launcher (Windows / macOS / Linux) for an already-installed Python environment.

    python run.py                    # uses config.ini next to this file
    python run.py my.ini --set space.max_heavy=9
    python run.py --selftest

On Linux servers prefer run.sh (it unpacks the bundled environment). Here any Python >= 3.10 with
rdkit, numpy and duckdb installed works (see ../environment.yml or requirements.txt).
"""
import os
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent


def main():
    os.environ.setdefault("PYTHONUNBUFFERED", "1")
    for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ.setdefault(v, "1")
    sys.path.insert(0, str(HERE))
    args = sys.argv[1:]
    if args and args[0] == "--selftest":
        import tempfile
        out = Path(tempfile.mkdtemp(prefix="ec_selftest_"))
        args = ["--config", str(HERE / "config.ini"), "--set", f"run.out={out}", "--set", "space.max_heavy=6",
                "--set", "space.allow_peroxide=no", "--set", "write.derive_max_heavy=6", "--set", "viewer.max_heavy=6",
                "--set", "run.workers=2", "--stages", "enumerate,derive,stats,viewer"]
        from ecspace import main as m
        m.main(args)
        import json
        got = [json.loads((out / "closed" / f"L{n}" / "DONE.json").read_text())["count"] for n in range(1, 7)]
        assert got == [2, 5, 17, 77, 379, 2270], got
        print("[selftest] 通过：", got, "输出目录", out)
        return
    cfg = str(HERE / "config.ini")
    if args and not args[0].startswith("--"):
        cfg = args.pop(0)
    from ecspace import main as m
    m.main(["--config", cfg] + args)


if __name__ == "__main__":
    main()
