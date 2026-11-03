"""Build a single self-contained HTML (data + smiles-drawer inlined) that opens from file:// offline.

python build_offline.py [data_dir] [out.html]

index.html reads <script id="data-nodes"> / <script id="data-reactions"> when they exist and only
falls back to fetch() otherwise, so nothing in the page logic changes.
"""
import sys
from pathlib import Path

here = Path(__file__).parent
data_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else here
out = Path(sys.argv[2]) if len(sys.argv) > 2 else here / "offline" / "CO2加氢_C1-C16反应网络.html"
src = (here / "index.html").read_text(encoding="utf-8")
lib = (here / "vendor" / "smiles-drawer.min.js").read_text(encoding="utf-8")


def tag(name):
    data = (data_dir / f"{name}.json").read_text(encoding="utf-8").replace("</", "<\\/")
    return f'<script id="data-{name}" type="application/json">{data}</script>'


src = "\n".join(l for l in src.splitlines() if "fonts.googleapis.com" not in l and "fonts.gstatic.com" not in l)
cdn = '<script src="https://cdn.jsdelivr.net/npm/smiles-drawer@2.4.1/dist/smiles-drawer.min.js"></script>'
assert cdn in src
src = src.replace(cdn, "<script>/* smiles-drawer 2.4.1, MIT licence (see smiles-drawer-LICENSE.md) */\n"
                  + lib.replace("</script", "<\\/script") + "\n</script>\n" + tag("nodes") + "\n" + tag("reactions"))
src = src.replace("<title>", "<!doctype html>\n<html lang=\"zh-CN\"><head><meta charset=\"utf-8\">\n"
                  "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1, viewport-fit=cover\">\n"
                  "<style>:root{color-scheme:light}body{margin:0}[hidden]{display:none!important}</style>\n<title>", 1)
src = src.replace("<div class=\"app\">", "</head><body>\n<div class=\"app\">", 1) + "\n</body></html>\n"
out.parent.mkdir(parents=True, exist_ok=True)
out.write_text(src, encoding="utf-8")
print(out, round(out.stat().st_size / 1e6, 2), "MB")
