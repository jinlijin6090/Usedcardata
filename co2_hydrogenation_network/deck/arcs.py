"""Arc diagram of olefin cracking / oligomerization fluxes (G1 calibrated network)."""
import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import Arc  # noqa: E402

d = json.load(open("fluxes.json"))["G1_GaZrOx_GaZSM5"]
ORANGE, BLUE, INK2, GRID = "#EB6834", "#2A78D6", "#52514E", "#E4E3DF"
fig, ax = plt.subplots(figsize=(11, 3.9), dpi=200)
fig.patch.set_facecolor("white")
ax.set_facecolor("white")
crack, olig = [], []
for r in d["reactions"]:
    st = r["stoich"]
    if r["family"] == "Z5_crack":
        src = [k for k, v in st.items() if v < 0][0]
        for k, v in st.items():
            if v > 0:
                crack.append((int(src[3:]), int(k[3:]), r["flux_rel"] * v))
    elif r["family"] == "Z4_olig":
        dst = [k for k, v in st.items() if v > 0][0]
        for k, v in st.items():
            if v < 0:
                olig.append((int(k[3:]), int(dst[3:]), r["flux_rel"] * (-v)))
fmax = max(f for *_, f in crack + olig)
for (a, b, f), above, col in [(x, True, ORANGE) for x in crack] + [(x, False, BLUE) for x in olig]:
    if f < 2e-4:
        continue
    cx, w = (a + b) / 2, abs(a - b)
    lw = 0.4 + 6 * (f / fmax) ** 0.7
    ax.add_patch(Arc((cx, 0), w, w * 0.9, theta1=0 if above else 180, theta2=180 if above else 360,
                     color=col, lw=lw, alpha=0.85))
ns = np.arange(2, 17)
ax.scatter(ns, np.zeros_like(ns), s=520, color="white", edgecolor=INK2, zorder=3, lw=1.2)
for n in ns:
    ax.text(n, 0, f"C{n}=", ha="center", va="center", fontsize=8, color="#0B0B0B", zorder=4)
ax.set_xlim(1.2, 16.8)
ax.set_ylim(-4.3, 5.7)
ax.axis("off")
fig.tight_layout()
fig.savefig("arcs.png", dpi=200, facecolor="white", bbox_inches="tight", pad_inches=0.05)
from PIL import Image
import json as _j
_j.dump(list(Image.open("arcs.png").size), open("arcs_size.json", "w"))
print("crack edges", len(crack), "olig edges", len(olig), "max flux", round(fmax, 4))
