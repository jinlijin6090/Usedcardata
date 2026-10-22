import matplotlib as mpl

# categorical order (validated reference palette) and an ordinal blue ramp for sizes
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
BLUE_RAMP = ["#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#104281", "#0d366b"]


def setup():
    mpl.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 10.5,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.edgecolor": "#8a8984",
        "axes.labelcolor": "#0b0b0b",
        "axes.titlesize": 11,
        "axes.titleweight": "bold",
        "axes.grid": True,
        "grid.color": "#e6e5e0",
        "grid.linewidth": 0.7,
        "xtick.color": "#52514e",
        "ytick.color": "#52514e",
        "lines.linewidth": 2,
        "figure.facecolor": "white",
        "savefig.bbox": "tight",
    })
