"""按化学类别（醚、烯醇、炔、乙烯酮、季碳、非甲基支链、卡宾等）给 species.csv 中的物种分类并统计。

    python tools/classify_species.py output_n12/species.csv classes.json
"""
import csv, sys, json, collections
from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")

P = {k: Chem.MolFromSmarts(v) for k, v in {
    "ether": "[#6]-[#8]-[#6]",
    "ester": "[#6]-[#8]-[#6]=[#8]",
    "acid": "[OX2H1]-[#6]=[#8]",
    "enol": "[#6]=[#6]-[OX2H1]",
    "enol_ether": "[#6]=[#6]-[#8]-[#6]",
    "alkyne": "[#6]#[#6]",
    "cumulene": "[#6]=[#6]=[#6]",
    "ketene": "[#6]=[#6]=[#8]",
    "quaternary_C": "[#6]([#6])([#6])([#6])[#6]",
    "gem_diol_like": "[#8]-[#6]-[#8]",
}.items()}

def carbon_long_branch(m):
    """碳骨架上是否存在 >=2 个碳的侧链（以最长碳链为主链）。"""
    cs = [a.GetIdx() for a in m.GetAtoms() if a.GetSymbol() == "C"]
    if len(cs) < 4:
        return False
    adj = {i: [n.GetIdx() for n in m.GetAtomWithIdx(i).GetNeighbors() if n.GetSymbol() == "C"] for i in cs}
    def bfs(s):
        dist, prev, q = {s: 0}, {s: None}, [s]
        for v in q:
            for w in adj[v]:
                if w not in dist:
                    dist[w] = dist[v] + 1; prev[w] = v; q.append(w)
        far = max(dist, key=dist.get)
        return far, dist, prev
    # 碳骨架可能被 O 隔开（醚），分别处理每个碳连通分量
    seen, long_branch = set(), False
    for c0 in cs:
        if c0 in seen:
            continue
        a, d, _ = bfs(c0); comp = set(d); seen |= comp
        b, d2, prev = bfs(a)
        path, v = set(), b
        while v is not None:
            path.add(v); v = prev[v]
        rest = comp - path
        # 剩余碳的连通分量大小 >= 2 即为非甲基支链
        done = set()
        for r in rest:
            if r in done:
                continue
            st, size = [r], 0; done.add(r)
            while st:
                v = st.pop(); size += 1
                for w in adj[v]:
                    if w in rest and w not in done:
                        done.add(w); st.append(w)
            if size >= 2:
                long_branch = True
    return long_branch

def classify(r):
    smi = r["smiles"]; heavy = int(r["heavy_atoms"]); nC = int(r["nC"]); nH = int(r["nH"]); nO = int(r["nO"])
    fv = int(r["free_valence"]); npi = int(r["n_pi"])
    m = Chem.MolFromSmiles(smi, sanitize=False)
    m.UpdatePropertyCache(strict=False)
    Chem.SanitizeMol(m, Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES, catchErrors=True)
    tags = set()
    has = lambda k: m.HasSubstructMatch(P[k])
    if nO == 2:
        os_ = [a for a in m.GetAtoms() if a.GetSymbol() == "O"]
        cn = [frozenset(n.GetIdx() for n in o.GetNeighbors() if n.GetSymbol() == "C") for o in os_]
        same_c = bool(cn[0] & cn[1])
        tags.add("two_O_same_C" if same_c else "two_O_different_C")
    if has("ester"): tags.add("ester")
    elif has("ether"): tags.add("ether")
    if has("acid"): tags.add("carboxylic_acid")
    if has("enol"): tags.add("enol")
    if has("enol_ether"): tags.add("enol_ether")
    if has("alkyne") and smi not in ("[C]#[C]",): tags.add("alkyne")
    if has("cumulene"): tags.add("cumulene/allene")
    if has("ketene"): tags.add("ketene")
    if has("quaternary_C"): tags.add("quaternary_C")
    if carbon_long_branch(m): tags.add("non_methyl_branch")
    if nC >= 2 and nH == 0: tags.add("no_H_carbon_cluster")
    if fv >= 2 and heavy >= 2 and smi != "[C]=O": tags.add("carbene_or_diradical(heavy>=2)")
    if fv == 0: tags.add("closed_shell")
    else: tags.add("open_shell")
    return tags

def main(path, out):
    tag_count = collections.Counter(); tag_small = collections.Counter(); ex = collections.defaultdict(list)
    n = 0; any_flag = 0
    for r in csv.DictReader(open(path)):
        n += 1
        tags = classify(r)
        for t in tags:
            tag_count[t] += 1
            if int(r["heavy_atoms"]) <= 5: tag_small[t] += 1
            if len(ex[t]) < 12: ex[t].append(r["smiles"])
    res = {"n_species": n, "counts": dict(tag_count), "counts_le5": dict(tag_small), "examples": ex}
    json.dump(res, open(out, "w"), indent=1, ensure_ascii=False)
    for t, c in tag_count.most_common():
        print(f"{t:34s} {c:8d} ({100*c/n:5.1f}%)  <=5: {tag_small[t]:5d}   e.g. {', '.join(ex[t][:5])}")

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
