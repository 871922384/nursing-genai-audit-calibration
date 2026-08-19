#!/usr/bin/env python3
"""Compute pilot interrater reliability without third-party dependencies."""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter, defaultdict
from pathlib import Path


RESULT_HEADER = (
    "variable_id",
    "domain",
    "label",
    "method",
    "critical",
    "n_total",
    "n_complete_pairs",
    "n_missing_pairs",
    "n_both_na",
    "n_applicability_disagreements",
    "percent_agreement",
    "cohen_kappa",
    "gwet_ac1",
    "weighted_kappa",
    "icc_a1",
    "prevalence_note",
    "target_status",
)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def paired_index(rows: list[dict[str, str]]) -> dict[tuple[str, str], dict[str, str]]:
    result: dict[tuple[str, str], dict[str, str]] = {}
    for row in rows:
        key = (row["pilot_id"], row["variable_id"])
        if key in result:
            raise ValueError(f"duplicate coding row: {key}")
        result[key] = row
    return result


def safe_ratio(numerator: float, denominator: float) -> float | None:
    return None if denominator == 0 else numerator / denominator


def cohen_kappa(pairs: list[tuple[str, str]]) -> float | None:
    if not pairs:
        return None
    categories = sorted({value for pair in pairs for value in pair})
    n = len(pairs)
    po = sum(a == b for a, b in pairs) / n
    left = Counter(a for a, _ in pairs)
    right = Counter(b for _, b in pairs)
    pe = sum(left[c] * right[c] for c in categories) / (n * n)
    return safe_ratio(po - pe, 1 - pe)


def gwet_ac1(pairs: list[tuple[str, str]]) -> float | None:
    if not pairs:
        return None
    categories = sorted({value for pair in pairs for value in pair})
    q = len(categories)
    if q < 2:
        return None
    n = len(pairs)
    po = sum(a == b for a, b in pairs) / n
    pooled = Counter(value for pair in pairs for value in pair)
    proportions = [pooled[c] / (2 * n) for c in categories]
    pe = sum(p * (1 - p) for p in proportions) / (q - 1)
    return safe_ratio(po - pe, 1 - pe)


def weighted_kappa(pairs: list[tuple[str, str]], order: list[str]) -> float | None:
    pairs = [(a, b) for a, b in pairs if a in order and b in order]
    if not pairs or len(order) < 2:
        return None
    rank = {value: index for index, value in enumerate(order)}
    denominator = len(order) - 1
    weight = lambda a, b: 1 - abs(rank[a] - rank[b]) / denominator
    n = len(pairs)
    observed = sum(weight(a, b) for a, b in pairs) / n
    left = Counter(a for a, _ in pairs)
    right = Counter(b for _, b in pairs)
    expected = sum(
        weight(a, b) * left[a] * right[b] / (n * n) for a in order for b in order
    )
    return safe_ratio(observed - expected, 1 - expected)


def icc_a1(pairs: list[tuple[str, str]]) -> float | None:
    try:
        values = [(float(a), float(b)) for a, b in pairs]
    except ValueError:
        return None
    n = len(values)
    k = 2
    if n < 2:
        return None
    grand = sum(sum(row) for row in values) / (n * k)
    row_means = [sum(row) / k for row in values]
    col_means = [sum(row[j] for row in values) / n for j in range(k)]
    ss_rows = k * sum((mean - grand) ** 2 for mean in row_means)
    ss_cols = n * sum((mean - grand) ** 2 for mean in col_means)
    ss_total = sum((value - grand) ** 2 for row in values for value in row)
    ss_error = ss_total - ss_rows - ss_cols
    ms_rows = ss_rows / (n - 1)
    ms_cols = ss_cols / (k - 1)
    ms_error = ss_error / ((n - 1) * (k - 1))
    denominator = ms_rows + (k - 1) * ms_error + k * (ms_cols - ms_error) / n
    return safe_ratio(ms_rows - ms_error, denominator)


def fmt(value: float | None) -> str:
    if value is None or not math.isfinite(value):
        return "NA"
    return f"{value:.4f}"


def target_status(
    *, critical: bool, agreement: float | None, kappa: float | None, ac1: float | None
) -> str:
    agreement_target = 0.85 if critical else 0.80
    coefficient_target = 0.80 if critical else 0.70
    if agreement is None:
        return "not_estimable"
    available = [value for value in (kappa, ac1) if value is not None and math.isfinite(value)]
    if not available:
        return "manual_review_required"
    return "pass" if agreement >= agreement_target and max(available) >= coefficient_target else "revise"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reviewer-a", required=True, type=Path)
    parser.add_argument("--reviewer-b", required=True, type=Path)
    parser.add_argument("--dictionary", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--adjudication-output", required=True, type=Path)
    args = parser.parse_args()

    dictionary = {row["variable_id"]: row for row in read_csv(args.dictionary)}
    index_a = paired_index(read_csv(args.reviewer_a))
    index_b = paired_index(read_csv(args.reviewer_b))
    if set(index_a) != set(index_b):
        raise ValueError("reviewer sheets do not contain identical pilot-variable keys")

    grouped: dict[str, list[tuple[str, str, str]]] = defaultdict(list)
    adjudication: list[dict[str, str]] = []
    for key in sorted(index_a):
        pilot_id, variable_id = key
        row_a, row_b = index_a[key], index_b[key]
        value_a, value_b = row_a["value"].strip(), row_b["value"].strip()
        grouped[variable_id].append((pilot_id, value_a, value_b))
        adjudication.append(
            {
                "pilot_id": pilot_id,
                "pmid": row_a["pmid"],
                "variable_id": variable_id,
                "domain": row_a["domain"],
                "label": row_a["label"],
                "reviewer_a_value": value_a,
                "reviewer_b_value": value_b,
                "agreement_before_adjudication": "" if not value_a or not value_b else ("Y" if value_a == value_b else "N"),
                "adjudicated_value": value_a if value_a and value_a == value_b else "",
                "adjudication_reason": "",
                "manual_change_required": "Y" if value_a != value_b or not value_a or not value_b else "N",
            }
        )

    results: list[dict[str, str]] = []
    for variable_id, meta in dictionary.items():
        rows = grouped.get(variable_id, [])
        complete = [(a, b) for _, a, b in rows if a and b]
        missing = sum(not a or not b for _, a, b in rows)
        both_na = sum(a == "NA" and b == "NA" for a, b in complete)
        applicability_disagreements = sum((a == "NA") != (b == "NA") for a, b in complete)
        nominal_pairs = [(a, b) for a, b in complete if not (a == "NA" and b == "NA")]
        agreement = safe_ratio(sum(a == b for a, b in nominal_pairs), len(nominal_pairs))
        method = meta["reliability_method"]
        kappa = ac1 = wkappa = icc = None
        if method == "Cohen_kappa_and_AC1":
            kappa = cohen_kappa(nominal_pairs)
            ac1 = gwet_ac1(nominal_pairs)
        elif method == "weighted_kappa_linear":
            order = [value for value in meta["allowed_values"].split("|") if value != "NA"]
            applicable_pairs = [(a, b) for a, b in nominal_pairs if a != "NA" and b != "NA"]
            wkappa = weighted_kappa(applicable_pairs, order)
            kappa = cohen_kappa(nominal_pairs)
            ac1 = gwet_ac1(nominal_pairs)
        elif method == "ICC_A1":
            applicable_pairs = [(a, b) for a, b in nominal_pairs if a != "NA" and b != "NA"]
            icc = icc_a1(applicable_pairs)

        pooled = Counter(value for pair in nominal_pairs for value in pair)
        top = pooled.most_common(1)[0] if pooled else ("", 0)
        top_share = safe_ratio(top[1], sum(pooled.values()))
        prevalence_note = ""
        if top_share is not None and top_share >= 0.90:
            prevalence_note = f"skewed:{top[0]}={top_share:.1%}; interpret AC1 with kappa"
        critical = meta["critical"] == "Y"
        coefficient_for_target = wkappa if method == "weighted_kappa_linear" else icc if method == "ICC_A1" else kappa
        alternate = ac1 if method == "Cohen_kappa_and_AC1" else coefficient_for_target
        results.append(
            {
                "variable_id": variable_id,
                "domain": meta["domain"],
                "label": meta["label"],
                "method": method,
                "critical": meta["critical"],
                "n_total": str(len(rows)),
                "n_complete_pairs": str(len(complete)),
                "n_missing_pairs": str(missing),
                "n_both_na": str(both_na),
                "n_applicability_disagreements": str(applicability_disagreements),
                "percent_agreement": fmt(agreement),
                "cohen_kappa": fmt(kappa),
                "gwet_ac1": fmt(ac1),
                "weighted_kappa": fmt(wkappa),
                "icc_a1": fmt(icc),
                "prevalence_note": prevalence_note,
                "target_status": target_status(
                    critical=critical,
                    agreement=agreement,
                    kappa=coefficient_for_target,
                    ac1=alternate,
                ) if method not in {"exact_match", "not_applicable"} else "not_scored",
            }
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=RESULT_HEADER)
        writer.writeheader()
        writer.writerows(results)
    with args.adjudication_output.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=tuple(adjudication[0]))
        writer.writeheader()
        writer.writerows(adjudication)

    summary = Counter(row["target_status"] for row in results)
    print(json.dumps({"variables": len(results), "target_status": summary}, ensure_ascii=False, default=dict))
    return 0 if all(row["n_missing_pairs"] == "0" for row in results) else 2


if __name__ == "__main__":
    raise SystemExit(main())
