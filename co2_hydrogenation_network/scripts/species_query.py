"""Query the species-level database (output of species_merge.py) with DuckDB, straight from the .tsv.gz files.

    python species_query.py <db_dir> --node "CC=C[CH]CC"            # one species: node, states, merged writings
    python species_query.py <db_dir> --reactions-of "CC(C)=C"        # every reaction of one species
    python species_query.py <db_dir> --old "CC[CH+]C"                # where an old (ion / R* / R.) key went
    python species_query.py <db_dir> "SELECT kind, count(*) FROM nodes GROUP BY 1 ORDER BY 2 DESC"
    python species_query.py <db_dir> --export out.parquet "SELECT * FROM nodes WHERE C = 8"

Views
  nodes(node_id, smiles, label, kind, formula, C, H, O, states, n_forms, skeleton)
  reactions(id, reactants, products, families, n_original, equation)  -- reactants / products are lists of node_id
  aliases(old_id, old_class, old_key, node_id, form_smiles)
SMILES given on the command line are canonicalised with RDKit first (any valid writing works; a "+"-marked
site key of the old database is looked up in aliases).
"""
import sys
from pathlib import Path

import duckdb


def connect(d):
    d = Path(d)
    con = duckdb.connect()
    rx = sorted(d.glob("reactions.part*.tsv.gz")) or [d / "reactions.tsv.gz"]
    opt = "delim='\\t', header=true, quote='', escape='', all_varchar=true"
    con.execute(f"""CREATE VIEW nodes AS SELECT node_id::BIGINT node_id, smiles, label, kind, formula, C::INT C, H::INT H,
                    O::INT O, states, n_forms::INT n_forms, skeleton FROM read_csv('{d / 'nodes.tsv.gz'}', {opt})""")
    con.execute(f"""CREATE VIEW reactions AS SELECT id::BIGINT id, string_split(reactants, ',')::BIGINT[] reactants,
                    string_split(products, ',')::BIGINT[] products, string_split(families, ';') families,
                    n_original::INT n_original, equation FROM read_csv({[str(p) for p in rx]!r}, {opt})""")
    con.execute(f"""CREATE VIEW aliases AS SELECT old_id::BIGINT old_id, old_class, old_key, TRY_CAST(node_id AS BIGINT) node_id, form_smiles
                    FROM read_csv('{d / 'aliases.tsv.gz'}', {opt})""")
    return con


def canon(s):
    try:
        from rdkit import Chem
        m = Chem.MolFromSmiles(s)
        return Chem.MolToSmiles(m) if m is not None else s
    except ImportError:
        return s


def show(cur):
    cols = [c[0] for c in cur.description]
    rows = cur.fetchall()
    w = [max([len(str(c))] + [len(str(r[i])) for r in rows]) for i, c in enumerate(cols)]
    print("  ".join(str(c).ljust(w[i]) for i, c in enumerate(cols)))
    print("  ".join("-" * x for x in w))
    for r in rows:
        print("  ".join(str(v).ljust(w[i]) for i, v in enumerate(r)))
    print(f"({len(rows)} 行)")


def node_ids(con, s):
    q = s.replace("'", "''")
    if "+" in s and "[" in s:          # old site key
        return [r[0] for r in con.execute(f"SELECT DISTINCT node_id FROM aliases WHERE old_key = '{q}'").fetchall()]
    c = canon(s).replace("'", "''")
    ids = [r[0] for r in con.execute(f"SELECT node_id FROM nodes WHERE smiles = '{c}' OR label = '{q}'").fetchall()]
    if not ids:                         # a resonance writing that is not the representative
        ids = [r[0] for r in con.execute(f"SELECT DISTINCT node_id FROM aliases WHERE form_smiles = '{c}'").fetchall()]
    return ids


def main():
    a = sys.argv[1:]
    if len(a) < 2:
        print(__doc__)
        return
    con = connect(a[0])
    if a[1] in ("--node", "--reactions-of", "--old"):
        ids = node_ids(con, a[2])
        if not ids:
            print("没有找到这个物种")
            return
        il = ",".join(map(str, ids))
        if a[1] == "--node":
            show(con.execute(f"SELECT * FROM nodes WHERE node_id IN ({il})"))
            show(con.execute(f"SELECT old_class, old_key, form_smiles FROM aliases WHERE node_id IN ({il}) ORDER BY 1, 2"))
        elif a[1] == "--old":
            show(con.execute(f"SELECT a.old_class, a.old_key, n.node_id, n.smiles, n.kind, n.states FROM aliases a "
                             f"JOIN nodes n USING (node_id) WHERE n.node_id IN ({il})"))
        else:
            show(con.execute(f"SELECT id, array_to_string(families, ';') families, n_original, equation FROM reactions "
                             f"WHERE list_has_any(reactants, [{il}]) OR list_has_any(products, [{il}]) ORDER BY id"))
    elif a[1] == "--export":
        con.execute(f"COPY ({a[3]}) TO '{a[2]}'")
        print("已写出", a[2])
    else:
        show(con.execute(a[1]))


if __name__ == "__main__":
    main()
