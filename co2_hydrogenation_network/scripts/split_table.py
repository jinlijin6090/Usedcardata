"""Split a large .tsv.gz into parts below a size limit (each part keeps the header line).

    python split_table.py /home/user/net_C16_full/reactions.tsv.gz ../data/net_C16_full/reactions 90
"""
import gzip
import sys
from pathlib import Path

src, prefix, limit_mb = sys.argv[1], sys.argv[2], float(sys.argv[3]) if len(sys.argv) > 3 else 90
limit = limit_mb * 1e6
with gzip.open(src, "rt") as fh:
    header = fh.readline()
    part, out, raw = 0, None, None
    for line in fh:
        if out is None or raw.tell() > limit:
            if out:
                out.close()
                raw.close()
            part += 1
            path = Path(f"{prefix}.part{part:02d}.tsv.gz")
            raw = open(path, "wb")
            out = gzip.GzipFile(fileobj=raw, mode="wb", compresslevel=6)
            out.write(header.encode())
        out.write(line.encode())
    out.close()
    raw.close()
print("parts", part)
