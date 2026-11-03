"""Collect the numbers the deck needs into deckdata.json."""
import json
import re

F = json.load(open("fluxes.json"))
g1 = F["G1_GaZrOx_GaZSM5"]


def nC(k):
    return 1 if k == "CH4" else int(re.sub(r"\D", "", k))


out = g1["outlet"]
hc = {k: v for k, v in out.items() if re.match(r"(ole|par|nap|aro)\d+$", k) or k == "CH4"}
tot = sum(v * nC(k) for k, v in hc.items())
share = {("par1" if k == "CH4" else k): v * nC(k) / tot for k, v in hc.items()}

edges = {}
def add(key, f):
    edges[key] = edges.get(key, 0) + f

for r in g1["reactions"]:
    st, fam, f = r["stoich"], r["family"], r["flux_rel"]
    neg = [k for k, v in st.items() if v < 0 and re.match(r"(ole|par|nap|aro)\d+$", k)]
    pos = [k for k, v in st.items() if v > 0 and re.match(r"(ole|par|nap|aro)\d+$", k)]
    if fam in ("Z3_meth", "Z8_ameth", "OX3_hyd", "Z6_cyc", "Z7b_dh", "Z7_HT") and neg and pos:
        add(f"{neg[0]}>{pos[0]}", f)
    if fam == "OX3_hyd" and neg:
        add(f"{neg[0]}>par{neg[0][3:]}", 0)  # ensure key
for k, v in g1["ht_partner_flux"].items():
    add(f"{k}>par{k[3:]}", v)          # paraffins made by hydride transfer

fam_order = [("MeOH_synthesis", "甲醇合成 OX1"), ("RWGS", "逆水煤气变换 OX2"), ("Z3_meth", "烯烃甲基化 Z3"),
             ("Z5_crack", "裂解 Z5"), ("Z8_ameth", "芳烃甲基化 Z8"), ("Z6_cyc", "环化 Z6"),
             ("Z7b_dh", "脱氢芳构化 Z7b"), ("OX3_hyd", "烯烃加氢 OX3"), ("Z9_dealk", "脱烷基 Z9"),
             ("Z7_HT", "氢转移芳构化 Z7"), ("Z4_olig", "齐聚 Z4"), ("Z10_ch4", "甲烷 Z10")]
cmp = {}
for case, d in F.items():
    vals = []
    for key, _ in fam_order:
        vals.append(d["oxide"].get(key, d["family_flux"].get(key, 0.0)))
    cmp[case] = vals

json.dump({"share": share, "edges": edges, "fam_labels": [x[1] for x in fam_order], "compare": cmp,
           "oxide": {c: d["oxide"] for c, d in F.items()},
           "family_flux": {c: d["family_flux"] for c, d in F.items()}},
          open("deckdata.json", "w"), indent=1, ensure_ascii=False)
print("edges", len(edges), "max", max(edges.values()))
