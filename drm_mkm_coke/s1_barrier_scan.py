"""Scan the three new barriers at the paper's base condition."""
from coke_tools import *
import itertools
pp = ProcessParams()
rows = []
# (a) Ea_cc x dH_cc grid, R19/R20 nominal
for Ea, dH in itertools.product(np.round(np.arange(1.0, 2.81, 0.1), 2), [-0.3, 0.0, 0.3]):
    cp = CokeParams(Ea_cc=Ea, dH_cc=dH)
    o = solve(load_model(cp), pp)
    rows.append(dict(scan="Ea_cc", Ea_cc=Ea, dH_cc=dH, Ea_ox=cp.Ea_ox, Ea_dep=cp.Ea_dep,
                     **{k: v for k, v in o.items() if not k.startswith("_")}))
# (b) Ea_dep scan at Ea_cc = 1.3 and 1.6
for Ea_cc in (1.3, 1.6):
    for Ed in np.round(np.arange(0.8, 2.41, 0.1), 2):
        cp = CokeParams(Ea_cc=Ea_cc, Ea_dep=Ed)
        o = solve(load_model(cp), pp)
        rows.append(dict(scan="Ea_dep", Ea_cc=Ea_cc, dH_cc=cp.dH_cc, Ea_ox=cp.Ea_ox, Ea_dep=Ed,
                         **{k: v for k, v in o.items() if not k.startswith("_")}))
# (c) Ea_ox scan
for Ea_cc in (1.3, 1.6):
    for Eo in np.round(np.arange(0.8, 2.01, 0.1), 2):
        cp = CokeParams(Ea_cc=Ea_cc, Ea_ox=Eo)
        o = solve(load_model(cp), pp)
        rows.append(dict(scan="Ea_ox", Ea_cc=Ea_cc, dH_cc=cp.dH_cc, Ea_ox=Eo, Ea_dep=cp.Ea_dep,
                         **{k: v for k, v in o.items() if not k.startswith("_")}))
df = pd.DataFrame(rows); df.to_csv(f"{RES}/s1_barrier_scan.csv", index=False)

a = df[df.scan == "Ea_cc"]
fig, ax = plt.subplots(1, 3, figsize=(13, 3.8))
for dH, g in a.groupby("dH_cc"):
    ax[0].semilogy(g.Ea_cc, g.S_coke.clip(1e-12), marker=".", label=f"ΔH_cc={dH:+.1f} eV")
    ax[1].semilogy(g.Ea_cc, g.r_coke.clip(1e-15), marker=".")
    ax[2].plot(g.Ea_cc, g.r_DRM, marker=".")
ax[0].set(xlabel="E$_a$(C*+C*→C$_2$*) [eV]", ylabel="coke selectivity  2r$_{20}$/r$_{CH_4}$", ylim=(1e-5, 2)); ax[0].legend(fontsize=8)
ax[1].set(xlabel="E$_a$(C*+C*→C$_2$*) [eV]", ylabel="coke rate [mol C g$^{-1}$ s$^{-1}$]", ylim=(1e-7, 1))
ax[2].set(xlabel="E$_a$(C*+C*→C$_2$*) [eV]", ylabel="DRM rate r$_{H_2}$ [mol g$^{-1}$ s$^{-1}$]")
for x in ax: x.axvline(1.5, color="gray", ls=":", lw=0.8)
fig.suptitle("Base condition (973 K, pCH4=pCO2=0.5 bar, r_m=4 nm): sensitivity to the C–C coupling step", fontsize=10)
fig.savefig(f"{FIG}/s1_Ea_cc_scan.png")

b = df[df.scan == "Ea_dep"]; c = df[df.scan == "Ea_ox"]
fig, ax = plt.subplots(1, 2, figsize=(9, 3.6))
for Ea_cc, g in b.groupby("Ea_cc"):
    ax[0].semilogy(g.Ea_dep, g.S_coke.clip(1e-12), marker=".", label=f"S_coke, Ea_cc={Ea_cc}")
    ax[0].semilogy(g.Ea_dep, g["th_C2*"].clip(1e-12), marker="x", ls="--", label=f"θ(C2*), Ea_cc={Ea_cc}")
ax[0].set(xlabel="E$_a$(C$_2$*→C$_{dep}$) [eV]", ylabel="coke selectivity / θ(C2*)", ylim=(1e-8, 2)); ax[0].legend(fontsize=7)
for Ea_cc, g in c.groupby("Ea_cc"):
    ax[1].semilogy(g.Ea_ox, g.S_coke.clip(1e-12), marker=".", label=f"Ea_cc={Ea_cc}")
ax[1].set(xlabel="E$_a$(C$_2$*+O*→CO*+C*) [eV]", ylabel="coke selectivity", ylim=(1e-8, 2)); ax[1].legend(fontsize=8)
fig.savefig(f"{FIG}/s1_Ea_dep_ox_scan.png")

pd.set_option("display.width", 200)
print(a[a.dH_cc == -0.3][["Ea_cc", "S_coke", "r_coke", "r_DRM", "th_C*", "th_C2*", "th_O*"]].to_string(index=False))
print(b[["Ea_cc", "Ea_dep", "S_coke", "th_C2*", "r_DRM"]].to_string(index=False))
print(c[["Ea_cc", "Ea_ox", "S_coke", "r_c2ox", "r_DRM"]].to_string(index=False))
