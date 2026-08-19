#!/usr/bin/env python3
"""Derive uncertainty intervals from the frozen RC5 paired coding values."""

from __future__ import annotations

import csv
import hashlib
import importlib.util
import math
from pathlib import Path
import random


TOPIC = Path(__file__).resolve().parents[1]
RELIABILITY = Path("data/rc5-reliability.csv")
REVIEWER_A = Path("data/not-distributed-reviewer-a.csv")
REVIEWER_B = Path("data/not-distributed-reviewer-b.csv")
DICTIONARY = Path("data/coding-dictionary.csv")


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def load_reliability_module(topic: Path):
    path = topic / "scripts/reliability_metrics.py"
    spec = importlib.util.spec_from_file_location("rc5_reliability", path)
    if spec is None or spec.loader is None:
        raise ValueError("cannot load reliability functions")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def wilson_interval(successes: int, total: int) -> tuple[float, float]:
    if total <= 0 or successes < 0 or successes > total:
        raise ValueError("invalid binomial counts")
    z = 1.959963984540054
    proportion = successes / total
    denominator = 1 + z * z / total
    center = (proportion + z * z / (2 * total)) / denominator
    half_width = (
        z
        * math.sqrt(
            proportion * (1 - proportion) / total + z * z / (4 * total * total)
        )
        / denominator
    )
    return center - half_width, center + half_width


def percentile(values: list[float], probability: float) -> float:
    if not values:
        raise ValueError("percentile requires at least one value")
    ordered = sorted(values)
    position = probability * (len(ordered) - 1)
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    fraction = position - lower
    return ordered[lower] * (1 - fraction) + ordered[upper] * fraction


def coefficient_value(module, method: str, pairs: list[tuple[str, str]], order: list[str]):
    if method == "Cohen_kappa_and_AC1":
        return module.gwet_ac1(pairs)
    if method == "weighted_kappa_linear":
        applicable = [(a, b) for a, b in pairs if a != "NA" and b != "NA"]
        return module.weighted_kappa(applicable, order)
    if method == "ICC_A1":
        applicable = [(a, b) for a, b in pairs if a != "NA" and b != "NA"]
        return module.icc_a1(applicable)
    return None


def bootstrap_interval(
    module,
    *,
    variable_id: str,
    method: str,
    pairs: list[tuple[str, str]],
    order: list[str],
    iterations: int,
) -> tuple[float | None, float | None, int]:
    if not pairs or iterations <= 0:
        return None, None, 0
    seed = int.from_bytes(hashlib.sha256(variable_id.encode()).digest()[:8], "big")
    generator = random.Random(seed)
    estimates: list[float] = []
    for _ in range(iterations):
        sample = [pairs[generator.randrange(len(pairs))] for _ in pairs]
        estimate = coefficient_value(module, method, sample, order)
        if estimate is not None and math.isfinite(estimate):
            estimates.append(estimate)
    if not estimates:
        return None, None, 0
    return percentile(estimates, 0.025), percentile(estimates, 0.975), len(estimates)


def fmt(value: float | None) -> str:
    return "NA" if value is None or not math.isfinite(value) else f"{value:.4f}"


def build_interval_rows(
    topic: Path = TOPIC, *, bootstrap_iterations: int = 10000
) -> list[dict[str, str]]:
    topic = topic.resolve()
    module = load_reliability_module(topic)
    reliability_path = topic / RELIABILITY
    dictionary_path = topic / DICTIONARY
    if not reliability_path.is_file():
        reliability_path = topic / "data/rc5-reliability.csv"
    if not dictionary_path.is_file():
        dictionary_path = topic / "data/coding-dictionary.csv"
    reliability = read_csv(reliability_path)
    dictionary = {row["variable_id"]: row for row in read_csv(dictionary_path)}
    reviewer_a_path = topic / REVIEWER_A
    reviewer_b_path = topic / REVIEWER_B
    if reviewer_a_path.is_file() and reviewer_b_path.is_file():
        reviewer_a = module.paired_index(read_csv(reviewer_a_path))
        reviewer_b = module.paired_index(read_csv(reviewer_b_path))
    else:
        combined = read_csv(topic / "data/paired-coding-values.csv")
        reviewer_a = module.paired_index(
            [
                {
                    "pilot_id": row["paper_id"],
                    "variable_id": row["variable_id"],
                    "value": row["coder_a_value"],
                }
                for row in combined
            ]
        )
        reviewer_b = module.paired_index(
            [
                {
                    "pilot_id": row["paper_id"],
                    "variable_id": row["variable_id"],
                    "value": row["coder_b_value"],
                }
                for row in combined
            ]
        )
    if set(reviewer_a) != set(reviewer_b):
        raise ValueError("reviewer sheets do not contain identical keys")

    grouped: dict[str, list[tuple[str, str]]] = {}
    for (pilot_id, variable_id), row_a in reviewer_a.items():
        del pilot_id
        row_b = reviewer_b[(row_a["pilot_id"], variable_id)]
        grouped.setdefault(variable_id, []).append(
            (row_a["value"].strip(), row_b["value"].strip())
        )

    enriched: list[dict[str, str]] = []
    for source in reliability:
        variable_id = source["variable_id"]
        complete = [(a, b) for a, b in grouped[variable_id] if a and b]
        pairs = [(a, b) for a, b in complete if not (a == "NA" and b == "NA")]
        agreements = sum(a == b for a, b in pairs)
        agreement_low, agreement_high = wilson_interval(agreements, len(pairs))
        method = source["method"]
        order = [
            value
            for value in dictionary[variable_id]["allowed_values"].split("|")
            if value != "NA"
        ]
        coefficient_name = {
            "Cohen_kappa_and_AC1": "gwet_ac1",
            "weighted_kappa_linear": "weighted_kappa",
            "ICC_A1": "icc_a1",
        }.get(method, "not_applicable")
        coefficient = coefficient_value(module, method, pairs, order)
        ci_low, ci_high, valid = bootstrap_interval(
            module,
            variable_id=variable_id,
            method=method,
            pairs=pairs,
            order=order,
            iterations=bootstrap_iterations,
        )
        row = dict(source)
        row.update(
            {
                "agreement_numerator": str(agreements),
                "agreement_denominator": str(len(pairs)),
                "agreement_ci_low": fmt(agreement_low),
                "agreement_ci_high": fmt(agreement_high),
                "decision_coefficient": coefficient_name,
                "decision_coefficient_value": fmt(coefficient),
                "decision_coefficient_ci_low": fmt(ci_low),
                "decision_coefficient_ci_high": fmt(ci_high),
                "bootstrap_valid_replicates": str(valid),
            }
        )
        enriched.append(row)
    return enriched


def strict_dual_coefficient_pass_count(rows: list[dict[str, str]]) -> int:
    count = 0
    for row in rows:
        agreement_target = 0.85 if row["critical"] == "Y" else 0.80
        coefficient_target = 0.80 if row["critical"] == "Y" else 0.70
        agreement = float(row["percent_agreement"])
        if row["method"] == "Cohen_kappa_and_AC1":
            values = (row["cohen_kappa"], row["gwet_ac1"])
            passes = all(
                value != "NA" and float(value) >= coefficient_target for value in values
            )
        elif row["method"] == "weighted_kappa_linear":
            value = row["weighted_kappa"]
            passes = value != "NA" and float(value) >= coefficient_target
        else:
            passes = False
        count += agreement >= agreement_target and passes
    return count


if __name__ == "__main__":
    rows = build_interval_rows()
    print(
        f"rows={len(rows)} strict_dual_coefficient_passes="
        f"{strict_dual_coefficient_pass_count(rows)}"
    )
