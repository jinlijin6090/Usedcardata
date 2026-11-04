"""Query the enumerated database with SQL (DuckDB reads the gzipped shards directly).

    python -m ecspace.query <out_dir> "SELECT type, count(*) FROM species GROUP BY 1 ORDER BY 2 DESC"
    python -m ecspace.query <out_dir> --reactions-of "O=C1OCCO1"
    python -m ecspace.query <out_dir> --export species.parquet "SELECT * FROM species WHERE heavy <= 8"

Views: species(key, type, charge, multiplicity, C, H, O, Li, heavy, parent)
       reactions(family, reactants, products, parent_heavy)
       closed(smiles, heavy)            -- every closed-shell molecule of every layer
"""
import sys
from pathlib import Path

import duckdb

SP_COLS = "{'key':'VARCHAR','type':'VARCHAR','charge':'INTEGER','multiplicity':'INTEGER','C':'INTEGER','H':'INTEGER','O':'INTEGER','Li':'INTEGER','heavy':'INTEGER','parent':'VARCHAR'}"
RX_COLS = "{'family':'VARCHAR','reactants':'VARCHAR','products':'VARCHAR','parent_heavy':'INTEGER'}"


def connect(out):
    out = Path(out)
    con = duckdb.connect()
    sp = [str(p) for p in sorted((out / "species").glob("L*/part-*.tsv.gz"))] + [str(out / "species" / "POOL.tsv")]
    rx = [str(p) for p in sorted((out / "reactions").glob("L*/part-*.tsv.gz"))]
    cl = sorted((out / "closed").glob("L*/s*.smi.gz"))
    if sp:
        con.execute(f"CREATE VIEW species AS SELECT * FROM read_csv({sp!r}, delim='\\t', header=false, columns={SP_COLS}, quote='', escape='', skip=0, null_padding=true) WHERE key <> 'key'")
    if rx:
        con.execute(f"CREATE VIEW reactions AS SELECT * FROM read_csv({rx!r}, delim='\\t', header=false, columns={RX_COLS}, quote='', escape='')")
    if cl:
        parts = " UNION ALL ".join(f"SELECT column0 AS smiles, {p.parent.name[1:]} AS heavy FROM read_csv('{p}', header=false, delim='\\t', quote='', escape='')" for p in cl)
        con.execute(f"CREATE VIEW closed AS {parts}")
    return con


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        return
    con = connect(sys.argv[1])
    if sys.argv[2] == "--reactions-of":
        k = sys.argv[3].replace("'", "''")
        q = f"SELECT family, reactants || ' ⇌ ' || products AS equation FROM reactions WHERE list_contains(string_split(reactants, ' + '), '{k}') OR list_contains(string_split(products, ' + '), '{k}') ORDER BY family"
    elif sys.argv[2] == "--export":
        dst, sql = sys.argv[3], sys.argv[4]
        con.execute(f"COPY ({sql}) TO '{dst}'")
        print("written", dst)
        return
    else:
        q = sys.argv[2]
    print(con.execute(q).df().to_string() if hasattr(con.execute(q), "df") and _has_pandas() else _fmt(con.execute(q)))


def _has_pandas():
    try:
        import pandas  # noqa: F401
        return True
    except ImportError:
        return False


def _fmt(cur):
    cols = [d[0] for d in cur.description]
    rows = cur.fetchall()
    w = [max(len(str(c)), *(len(str(r[i])) for r in rows)) if rows else len(str(c)) for i, c in enumerate(cols)]
    lines = ["  ".join(str(c).ljust(w[i]) for i, c in enumerate(cols)), "  ".join("-" * x for x in w)]
    lines += ["  ".join(str(v).ljust(w[i]) for i, v in enumerate(r)) for r in rows]
    return "\n".join(lines)


if __name__ == "__main__":
    main()
