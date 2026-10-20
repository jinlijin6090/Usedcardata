"""比较两次运行的网络是否相同（按物种规范键，与 id 编号无关）。

    python tools/compare_runs.py run_a/ run_b/
"""
import csv, gzip, os, sys, hashlib

def load(d):
    keys = {r['species_id']: r['connectivity_key'] for r in csv.DictReader(open(os.path.join(d, 'species.csv')))}
    p = os.path.join(d, 'reactions.csv.gz')
    fh = gzip.open(p, 'rt', newline='') if os.path.exists(p) else open(os.path.join(d, 'reactions.csv'), newline='')
    rx, n = set(), 0
    for r in csv.DictReader(fh):
        t = (r['type'], tuple(sorted(keys[i] for i in r['reactant_ids'].split())),
             tuple(sorted(keys[i] for i in r['product_ids'].split())))
        rx.add(hashlib.blake2b(repr(t).encode(), digest_size=8).digest()); n += 1
    return set(keys.values()), rx, n

a, b = load(sys.argv[1]), load(sys.argv[2])
print('species:', len(a[0]), len(b[0]), 'equal' if a[0] == b[0] else 'DIFFERENT')
print('reactions (rows):', a[2], b[2], ' unique:', len(a[1]), len(b[1]), 'equal' if a[1] == b[1] else 'DIFFERENT',
      ' only-in-first:', len(a[1] - b[1]), ' only-in-second:', len(b[1] - a[1]))
