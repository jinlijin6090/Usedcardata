"""Turn a LOD tile directory into a double-click offline package (works from file://).

    python lod_build.py ../data/net_C16_full pkg/data --js   # tiles + .js twins
    python lod_package.py pkg                                 # index page, library, readme
"""
import sys
from pathlib import Path

here = Path(__file__).parent
pkg = Path(sys.argv[1])
data = pkg / "data"
src = (here / "lod.html").read_text(encoding="utf-8")
src = "\n".join(l for l in src.splitlines() if "fonts.googleapis.com" not in l and "fonts.gstatic.com" not in l)
cdn = "https://cdn.jsdelivr.net/npm/smiles-drawer@2.4.1/dist/smiles-drawer.min.js"
assert cdn in src
src = src.replace(cdn, "smiles-drawer.min.js")
src = src.replace("<title>", '<!doctype html>\n<html lang="zh-CN"><head><meta charset="utf-8">\n'
                  '<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">\n'
                  "<style>:root{color-scheme:light}body{margin:0}[hidden]{display:none!important}</style>\n<title>", 1)
src = src.replace('<div class="app">', '</head><body>\n<div class="app">', 1) + "\n</body></html>\n"
(data / "index.html").write_text(src, encoding="utf-8")
(data / "smiles-drawer.min.js").write_bytes((here / "vendor" / "smiles-drawer.min.js").read_bytes())
(data / "smiles-drawer-LICENSE.md").write_bytes((here / "vendor" / "smiles-drawer-LICENSE.md").read_bytes())
# remove the online-only .bin twins from the package to halve its size
for kind in ("l", "x"):
    for f in (data / kind).glob("*.bin"):
        f.unlink()
(data / "meta.json").unlink(missing_ok=True)
launcher = '<!doctype html><meta charset="utf-8"><meta http-equiv="refresh" content="0; url=data/index.html"><a href="data/index.html">打开反应网络</a>'
(pkg / "打开反应网络.html").write_text(launcher, encoding="utf-8")
print("package ready:", pkg)
