"""Parse the supplementary-material tables (dumped to data/si_tables_raw.json)
into tidy CSVs used by the TEA model.

Outputs
  data/equipment_costs.csv : route, item, EU_EUR, US_USD, ME_SAR
  data/varopex_si.csv      : route, region, type, category, value, currency
"""
import csv
import json
import re
from pathlib import Path

HERE = Path(__file__).parent
RAW = json.loads((HERE / "data/si_tables_raw.json").read_text())

EQUIP_TABLES = {29: "FTS", 30: "MtJ-1", 31: "MtJ-2", 32: "EtJ"}
VAROPEX_TABLES = {
    37: ("FTS", "EU"), 38: ("MtJ-1", "EU"), 39: ("MtJ-2", "EU"), 40: ("EtJ", "EU"),
    41: ("FTS", "US"), 42: ("MtJ-1", "US"), 43: ("MtJ-2", "US"), 44: ("EtJ", "US"),
    45: ("FTS", "ME"), 46: ("MtJ-1", "ME"), 47: ("MtJ-2", "ME"), 48: ("EtJ", "ME"),
}


def num(s: str) -> float:
    s = re.sub(r"[^\d.\-]", "", s.replace(",", ""))
    return float(s) if s not in ("", "-", ".") else 0.0


def main() -> None:
    with open(HERE / "data/equipment_costs.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["route", "item", "EU_EUR", "US_USD", "ME_SAR"])
        for idx, route in EQUIP_TABLES.items():
            for row in RAW[idx][2:]:
                w.writerow([route, row[0], num(row[1]), num(row[2]), num(row[3])])

    with open(HERE / "data/varopex_si.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["route", "region", "type", "category", "value"])
        for idx, (route, region) in VAROPEX_TABLES.items():
            for row in RAW[idx][2:]:
                typ, cat, val = row[0], row[1], row[2]
                if typ.lower().startswith("variable"):
                    w.writerow([route, region, "TOTAL", "variable Opex", num(row[-1])])
                    continue
                if len(set(row)) < 3:
                    continue
                w.writerow([route, region, typ, cat, num(val)])


if __name__ == "__main__":
    main()
