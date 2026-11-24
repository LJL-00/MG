#!/usr/bin/env python3
"""Summarize MIC distributions and apply the stated interpretation thresholds.

Only the Python standard library is required. Values remain concentration bounds.
MIC50/90 use nearest ranks ceil(p*n), with the AB/FC lower-bound pooling specified
in the manuscript. Caspofungin categories are secondary reference categories.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import re
from collections import Counter
from dataclasses import dataclass
from decimal import Decimal
from pathlib import Path

DRUGS = ("AB", "AND", "MF", "CAS", "FC", "PZ", "VOR", "IZ", "FZ")
NAMES = dict(zip(DRUGS, ("Amphotericin B", "Anidulafungin", "Micafungin", "Caspofungin",
    "Flucytosine", "Posaconazole", "Voriconazole", "Itraconazole", "Fluconazole")))
ECVS = {"FZ": Decimal("8"), "IZ": Decimal("2"), "PZ": Decimal("0.5")}


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def write_csv(path, fields, rows):
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def number(v):
    return format(v.normalize(), "f")


@dataclass(frozen=True)
class MIC:
    operator: str
    bound: Decimal

    @classmethod
    def parse(cls, text):
        match = re.fullmatch(r"\s*(<=|>=|[<>≤≥=]?)\s*(\d+(?:\.\d+)?)\s*", text)
        if not match:
            raise ValueError(f"Invalid MIC: {text!r}")
        op = {"<=": "≤", ">=": "≥", "": "="}.get(match[1], match[1])
        bound = Decimal(match[2])
        if not bound.is_finite() or bound <= 0:
            raise ValueError(f"MIC must be positive: {text!r}")
        return cls(op, bound)

    def text(self):
        return ("" if self.operator == "=" else self.operator) + number(self.bound)

    def sort_key(self):
        return self.bound, {"<": 0, "≤": 1, "=": 2, "≥": 3, ">": 4}[self.operator]

    def ecv_category(self, cutoff):
        if self.operator in ("<", "≤", "=") and self.bound <= cutoff:
            return "WT"
        if self.operator in ("=", "≥") and self.bound > cutoff:
            return "NWT"
        if self.operator == ">" and self.bound >= cutoff:
            return "NWT"
        raise ValueError(f"MIC {self.text()} crosses ECV {cutoff}; no exact category can be assigned")

    def breakpoint_category(self):
        if self.operator in ("<", "≤", "=") and self.bound <= 2:
            return "S"
        if self.operator in ("=", "≥", ">") and self.bound >= 8:
            return "R"
        if self.operator == "=" and self.bound == 4:
            return "I"
        raise ValueError(f"MIC {self.text()} does not identify a single S/I/R category")


def summary_mic(drug, raw):
    mic = MIC.parse(raw)
    # Pool the bottom-dilution readings used for the MIC summary.
    limit = {"AB": Decimal("0.12"), "FC": Decimal("0.06")}.get(drug)
    if limit is not None and mic.bound <= limit and mic.operator in ("<", "≤", "="):
        return MIC("≤", limit)
    return mic


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", type=Path, required=True, help="MIC CSV")
    p.add_argument("--qc", type=Path, help="Optional QC CSV")
    p.add_argument("--output", dest="output_dir", type=Path, required=True)
    args = p.parse_args()
    output = args.output_dir
    output.mkdir(parents=True, exist_ok=True)
    rows = read_csv(args.input)
    if not rows:
        raise ValueError("The MIC CSV is empty")
    if len({r["isolate"] for r in rows}) != len(rows):
        raise ValueError("Duplicate isolate in MIC inputs")
    for r in rows:
        for drug in DRUGS:
            MIC.parse(r[drug])
    main_rows = [r for r in rows if r["analysis_set"] == "Main"]
    separate = [r for r in rows if r["analysis_set"] == "ANI check only"]
    if len(main_rows) + len(separate) != len(rows):
        raise ValueError("Unknown analysis set")
    if not main_rows:
        raise ValueError("At least one isolate must have analysis_set Main")
    classified = []
    for r in main_rows:
        rec = {k: r.get(k, "") for k in ("isolate", "analysis_set", "clade", *DRUGS)}
        for drug, cutoff in ECVS.items():
            rec[drug + "_category"] = MIC.parse(r[drug]).ecv_category(cutoff)
        for drug in ("AND", "MF", "CAS"):
            rec[drug + "_category"] = MIC.parse(r[drug]).breakpoint_category()
        rec["CAS_category_role"] = "secondary reference only"
        rec["joint_azole_nwt"] = "Yes" if all(rec[d + "_category"] == "NWT" for d in ECVS) else "No"
        classified.append(rec)
    write_csv(output / "isolate_classifications.csv", list(classified[0]), classified)
    write_csv(output / "separately_analysed_mics.csv", list(rows[0]), separate)
    summaries, raw_distribution, pooled_distribution = [], [], []
    for drug in DRUGS:
        values = sorted((summary_mic(drug, r[drug]) for r in main_rows), key=lambda x: x.sort_key())
        n = len(values)
        rank50, rank90 = math.ceil(0.5 * n), math.ceil(0.9 * n)
        counts = Counter(r.get(drug + "_category", "Not assigned") for r in classified)
        role = "ECV" if drug in ECVS else "clinical breakpoint" if drug in ("AND", "MF") else "secondary reference only" if drug == "CAS" else "MIC distribution only"
        rec = {"drug_code": drug, "drug": NAMES[drug], "n": n, "minimum": values[0].text(),
            "maximum": values[-1].text(), "MIC50": values[rank50-1].text(), "MIC90": values[rank90-1].text(),
            "MIC50_rank": rank50, "MIC90_rank": rank90, "classification_role": role,
            "ECV_ug_per_mL": number(ECVS[drug]) if drug in ECVS else "",
            **{k: counts.get(k, 0) if (drug in ECVS and k in ("WT", "NWT")) or (drug in ("AND", "MF", "CAS") and k in ("S", "I", "R")) else ""
               for k in ("WT", "NWT", "S", "I", "R")}}
        summaries.append(rec)
        for text, count in sorted(Counter(r[drug] for r in main_rows).items(), key=lambda x: MIC.parse(x[0]).sort_key()):
            raw_distribution.append({"drug_code": drug, "recorded_MIC": text, "n": count})
        for text, count in sorted(Counter(x.text() for x in values).items(), key=lambda x: MIC.parse(x[0]).sort_key()):
            pooled_distribution.append({"drug_code": drug, "summary_MIC": text, "n": count})
    write_csv(output / "mic_summary.csv", list(summaries[0]), summaries)
    write_csv(output / "mic_distribution_raw.csv", list(raw_distribution[0]), raw_distribution)
    write_csv(output / "mic_distribution_summary.csv", list(pooled_distribution[0]), pooled_distribution)
    if args.qc:
        qc = read_csv(args.qc)
        if not qc:
            raise ValueError("The QC CSV is empty")
        for r in qc:
            value = MIC.parse(r["final_mic"])
            if value.operator != "=":
                raise ValueError("A bounded QC MIC needs interval-specific acceptance handling")
            lower, upper = Decimal(r["lower_limit"]), Decimal(r["upper_limit"])
            if lower > upper:
                raise ValueError("QC lower limit exceeds upper limit")
            r["computed_in_range"] = "Yes" if lower <= value.bound <= upper else "No"
        write_csv(output / "qc_verification.csv", list(qc[0]), qc)
    joint = [r["isolate"] for r in classified if r["joint_azole_nwt"] == "Yes"]
    summary = {"units": "ug/mL", "study_isolates": len(rows), "main_isolates": len(main_rows),
        "separately_analysed_isolates": len(separate), "MIC_quantile_method": "nearest rank ceil(p*n)",
        "summary_pooling": {"AB": "<0.12, ≤0.12 and 0.12 pooled as ≤0.12", "FC": "<0.06 and ≤0.06 pooled as ≤0.06"},
        "classification": {"azole_ECVs": {k: float(v) for k,v in ECVS.items()},
            "echinocandins": "S <=2, I 4, R >=8", "CAS_role": "secondary reference only; excluded from main phenotype classification",
            "unclassified_drugs": ["AB", "FC", "VOR"]},
        "joint_azole_NWT_isolates": joint, "drug_summaries": summaries}
    (output / "ast_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"main_isolates": len(main_rows), "separate_isolates": len(separate), "joint_NWT": len(joint)}))



if __name__ == "__main__":
    main()
