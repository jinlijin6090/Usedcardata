"""Degree of rate control for DRM rate, coke rate and coke selectivity."""
from coke_tools import *
rows = []
for Ea in [1.5, 2.2]:
    rc = load_model(CokeParams(Ea_cc=Ea))
    for cond_name, pp in [("base", ProcessParams()), ("pCH4=1.5", pp_with(pCH4=1.5)), ("rm=10nm", pp_with(rm=10.0))]:
        for key in ["r_DRM", "r_coke", "S_coke"]:
            X = drc(rc, pp, key)
            for rxn, v in X.items():
                rows.append(dict(Ea_cc=Ea, condition=cond_name, objective=key, reaction=rxn, X=v))
df = pd.DataFrame(rows); df.to_csv(f"{RES}/s3_drc.csv", index=False)
piv = df.pivot_table(index=["Ea_cc", "condition", "reaction"], columns="objective", values="X").reset_index()
piv["rn"] = piv.reaction.str[1:].astype(int); piv = piv.sort_values(["Ea_cc", "condition", "rn"])
pd.set_option("display.width", 200)
print(piv.to_string(index=False, float_format=lambda v: f"{v:+.3f}"))
for Ea in [1.5, 2.2]:
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.8), sharey=True)
    for a, cond in zip(ax, ["base", "pCH4=1.5", "rm=10nm"]):
        g = piv[(piv.Ea_cc == Ea) & (piv.condition == cond)]
        y = np.arange(len(g)); w = 0.27
        a.barh(y - w, g.r_DRM, w, label="X(DRM rate)")
        a.barh(y, g.r_coke, w, label="X(coke rate)")
        a.barh(y + w, g.S_coke, w, label="X(coke selectivity)")
        a.set_yticks(y); a.set_yticklabels(g.reaction); a.axvline(0, color="k", lw=0.6); a.set_title(cond)
        a.set_xlabel("degree of rate control")
    ax[0].legend(fontsize=8); fig.suptitle(f"DRC, Ea_cc = {Ea} eV", fontsize=10); fig.tight_layout()
    fig.savefig(f"{FIG}/s3_drc_Ea{Ea}.png"); plt.close(fig)
