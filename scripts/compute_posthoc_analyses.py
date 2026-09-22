#!/usr/bin/env python3
"""Compute post hoc statistical analyses for RINAH transfer revision.

Analyses:
- A2: Cluster bootstrap 95% CIs for pooled Stage 1 and Stage 2 agreement (by paper, 10,000 resamples)
- A3: Stage 1 chance-corrected agreement (AC1 and Cohen's kappa for NA vs Applicable) across 26 items and pooled
- A4: Locator-based mechanism classification for the 75 calibration disagreements
- A5: Sensitivity analysis merging N and NR categories into 'not reported as performed'
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
import random
import re

REPO_ROOT = Path(__file__).resolve().parents[1]
PAIRED_CSV = REPO_ROOT / "data/paired-coding-values.csv"
LOCATORS_CSV = REPO_ROOT / "data/disagreement-locators.csv"
ITEM_MAP_CSV = REPO_ROOT / "data/item-id-map.csv"


def normalize_section_name(s: str) -> str:
    if not s:
        return ""
    s = s.replace("#", "").strip().lower()
    s = re.sub(r"\(.*?\)", "", s).strip()
    return s


def extract_candidate_sections(loc: str) -> list[str]:
    if not loc or "chunk" in loc.lower():
        return []
    clauses = loc.split(";")
    secs = []
    for cl in clauses:
        parts = re.split(r"[|>]", cl)
        if len(parts) >= 2:
            sec = normalize_section_name(parts[1])
            if sec and not re.match(r"^l\d+", sec):
                secs.append(sec)
        else:
            c_parts = cl.split(",")
            if len(c_parts) >= 2:
                sec = normalize_section_name(c_parts[1])
                if sec and not re.match(r"^l\d+", sec):
                    secs.append(sec)
            else:
                sec = normalize_section_name(cl)
                if sec and not re.match(r"^l\d+", sec):
                    secs.append(sec)
    return secs


def sections_match(s_a: str, s_b: str) -> bool:
    if s_a == s_b:
        return True
    if "limitations" in s_a and "limitations" in s_b:
        return True
    if "outcome measures" in s_a and "outcome measures" in s_b:
        return True
    if "chatbot platform" in s_a and "chatbot platform" in s_b:
        return True
    return False


def classify_locator_pair(loc_a: str, loc_b: str) -> str:
    """Deterministically categorize evidence locator discrepancies into three classes:
    1. Full-text absence: one rater identified an explicit passage while the other found none across the text (chunk coverage)
    2. Same section: both raters located candidate evidence within the same manuscript section
    3. Different sections: raters located candidate evidence across different manuscript sections
    """
    if not loc_a or not loc_b or "chunk" in loc_a.lower() or "chunk" in loc_b.lower():
        return "Full-text absence"

    secs_a = extract_candidate_sections(loc_a)
    secs_b = extract_candidate_sections(loc_b)

    for sa in secs_a:
        for sb in secs_b:
            if sections_match(sa, sb):
                return "Same section"
    return "Different sections"


def safe_ratio(num: float, den: float) -> float | None:
    return None if den == 0 else num / den


def cohen_kappa(pairs: list[tuple[str, str]]) -> float | None:
    if not pairs:
        return None
    categories = sorted({val for pair in pairs for val in pair})
    n = len(pairs)
    po = sum(a == b for a, b in pairs) / n
    left = Counter(a for a, _ in pairs)
    right = Counter(b for _, b in pairs)
    pe = sum(left[c] * right[c] for c in categories) / (n * n)
    return safe_ratio(po - pe, 1 - pe)


def gwet_ac1(pairs: list[tuple[str, str]]) -> float | None:
    if not pairs:
        return None
    categories = sorted({val for pair in pairs for val in pair})
    q = len(categories)
    if q < 2:
        return 1.0 if all(a == b for a, b in pairs) else None
    n = len(pairs)
    po = sum(a == b for a, b in pairs) / n
    pooled = Counter(val for pair in pairs for val in pair)
    proportions = [pooled[c] / (2 * n) for c in categories]
    pe = sum(p * (1 - p) for p in proportions) / (q - 1)
    return safe_ratio(po - pe, 1 - pe)


def wilson_ci(k: int, n: int, z: float = 1.959963984540054) -> tuple[float, float]:
    if n == 0:
        return (0.0, 0.0)
    p = k / n
    denom = 1 + z * z / n
    center = (p + z * z / (2 * n)) / denom
    hw = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
    return (max(0.0, center - hw), min(1.0, center + hw))


def percentile(vals: list[float], p: float) -> float:
    if not vals:
        raise ValueError("empty list")
    s = sorted(vals)
    pos = p * (len(s) - 1)
    lo = math.floor(pos)
    hi = math.ceil(pos)
    if lo == hi:
        return s[lo]
    return s[lo] * (1 - (pos - lo)) + s[hi] * (pos - lo)


def main():
    # 1. Load paired coding values
    with PAIRED_CSV.open(encoding="utf-8") as f:
        rows = list(csv.DictReader(f))

    # Group by paper_id
    papers = sorted({r["paper_id"] for r in rows})
    variables = sorted({r["variable_id"] for r in rows})

    print(f"Loaded {len(rows)} paired coding rows across {len(papers)} papers and {len(variables)} variables.")

    paper_data: dict[str, list[dict[str, str]]] = defaultdict(list)
    for r in rows:
        paper_data[r["paper_id"]].append(r)

    # Verify locked invariants
    n_total = len(rows)
    both_na = sum(1 for r in rows if r["coder_a_value"] == "NA" and r["coder_b_value"] == "NA")
    ad_count = sum(1 for r in rows if (r["coder_a_value"] == "NA") ^ (r["coder_b_value"] == "NA"))
    ja_count = sum(1 for r in rows if r["coder_a_value"] != "NA" and r["coder_b_value"] != "NA")
    s1_agree = sum(1 for r in rows if (r["coder_a_value"] == "NA") == (r["coder_b_value"] == "NA"))
    s2_agree = sum(1 for r in rows if r["coder_a_value"] != "NA" and r["coder_b_value"] != "NA" and r["coder_a_value"] == r["coder_b_value"])
    vd_count = ja_count - s2_agree

    assert n_total == 312, f"Expected 312, got {n_total}"
    assert both_na == 102, f"Expected 102, got {both_na}"
    assert ad_count == 47, f"Expected 47, got {ad_count}"
    assert ja_count == 163, f"Expected 163, got {ja_count}"
    assert s1_agree == 265, f"Expected 265, got {s1_agree}"
    assert s2_agree == 135, f"Expected 135, got {s2_agree}"
    assert vd_count == 28, f"Expected 28, got {vd_count}"
    print("Locked invariants strictly verified!")

    # -------------------------------------------------------------
    # A2: Cluster Bootstrap 95% CIs (by paper, 10,000 resamples)
    # -------------------------------------------------------------
    rng = random.Random(20260922)
    boot_s1_rates = []
    boot_s2_rates = []

    for _ in range(10000):
        sample_papers = [rng.choice(papers) for _ in range(len(papers))]
        s1_num = 0
        s1_den = 0
        s2_num = 0
        s2_den = 0
        for pid in sample_papers:
            for r in paper_data[pid]:
                s1_den += 1
                a_na = (r["coder_a_value"] == "NA")
                b_na = (r["coder_b_value"] == "NA")
                if a_na == b_na:
                    s1_num += 1
                if not a_na and not b_na:
                    s2_den += 1
                    if r["coder_a_value"] == r["coder_b_value"]:
                        s2_num += 1
        boot_s1_rates.append(s1_num / s1_den)
        if s2_den > 0:
            boot_s2_rates.append(s2_num / s2_den)

    s1_point = s1_agree / n_total
    s2_point = s2_agree / ja_count
    s1_ci_lo = percentile(boot_s1_rates, 0.025)
    s1_ci_hi = percentile(boot_s1_rates, 0.975)
    s2_ci_lo = percentile(boot_s2_rates, 0.025)
    s2_ci_hi = percentile(boot_s2_rates, 0.975)

    print("\n=== A2: CLUSTER BOOTSTRAP 95% CIs (by paper, B=10,000) ===")
    print(f"Stage 1 Applicability Agreement: {s1_point*100:.1f}% ({s1_agree}/{n_total}), 95% Cluster CI: [{s1_ci_lo*100:.1f}%, {s1_ci_hi*100:.1f}%]")
    print(f"Stage 2 Value Agreement:         {s2_point*100:.1f}% ({s2_agree}/{ja_count}), 95% Cluster CI: [{s2_ci_lo*100:.1f}%, {s2_ci_hi*100:.1f}%]")

    # -------------------------------------------------------------
    # A3: Stage 1 Chance-Corrected Agreement across 26 items and pooled
    # -------------------------------------------------------------
    # For Stage 1: each rater assigns "NA" or "APP"
    var_rows: dict[str, list[dict[str, str]]] = defaultdict(list)
    for r in rows:
        var_rows[r["variable_id"]].append(r)

    print("\n=== A3: STAGE 1 CHANCE-CORRECTED AGREEMENT (26 Items) ===")
    s1_item_results = {}
    pooled_s1_pairs = []
    for var in sorted(var_rows.keys()):
        v_list = var_rows[var]
        pairs = [("NA" if r["coder_a_value"] == "NA" else "APP",
                  "NA" if r["coder_b_value"] == "NA" else "APP") for r in v_list]
        pooled_s1_pairs.extend(pairs)
        n = len(pairs)
        agree_cnt = sum(a == b for a, b in pairs)
        raw_agree = agree_cnt / n
        k = cohen_kappa(pairs)
        ac1 = gwet_ac1(pairs)
        s1_item_results[var] = {
            "n": n,
            "agree": agree_cnt,
            "raw_pct": raw_agree * 100,
            "kappa": k,
            "ac1": ac1,
        }

    pooled_s1_raw = sum(a == b for a, b in pooled_s1_pairs) / len(pooled_s1_pairs)
    pooled_s1_k = cohen_kappa(pooled_s1_pairs)
    pooled_s1_ac1 = gwet_ac1(pooled_s1_pairs)
    print(f"Pooled Stage 1 (312 pairs): Raw Agree = {pooled_s1_raw*100:.1f}%, Cohen's kappa = {pooled_s1_k:.3f}, Gwet's AC1 = {pooled_s1_ac1:.3f}")

    # -------------------------------------------------------------
    # A4: Deterministic Classification of the 75 calibration disagreements
    # -------------------------------------------------------------
    loc_classes = Counter()
    ad_loc_classes = Counter()
    vd_loc_classes = Counter()
    with LOCATORS_CSV.open(encoding="utf-8-sig") as f:
        loc_rows = list(csv.DictReader(f))

    assert len(loc_rows) == 75, f"Expected 75 disagreement rows, got {len(loc_rows)}"

    for r in loc_rows:
        # Deterministically compute classification from raw locators
        computed_pat = classify_locator_pair(r["coder_a_locator"], r["coder_b_locator"])
        archived_pat = r["locator_pattern"]
        assert computed_pat == archived_pat, (
            f"Deterministic classification mismatch for {r['paper_id']} {r['variable_id']}: "
            f"computed '{computed_pat}' vs archived '{archived_pat}'"
        )
        is_ad = r["disagreement_type"] == "AD"
        loc_classes[computed_pat] += 1
        if is_ad:
            ad_loc_classes[computed_pat] += 1
        else:
            vd_loc_classes[computed_pat] += 1

    print("\n=== A4: LOCATOR MECHANISM CLASSIFICATION (75 Disagreements) ===")
    print(f"Total Disagreements: 75 (47 AD + 28 VD)")
    for cat, cnt in loc_classes.most_common():
        print(f"  {cat}: {cnt} total (AD: {ad_loc_classes[cat]}, VD: {vd_loc_classes[cat]})")

    # -------------------------------------------------------------
    # A5: Sensitivity Analysis merging N and NR
    # -------------------------------------------------------------
    # Merge N and NR into 'NOT_PERF' among jointly applicable pairs
    ja_pairs_original = []
    ja_pairs_merged = []
    for r in rows:
        a = r["coder_a_value"]
        b = r["coder_b_value"]
        if a != "NA" and b != "NA":
            ja_pairs_original.append((a, b))
            a_m = "NOT_PERF" if a in ("N", "NR") else a
            b_m = "NOT_PERF" if b in ("N", "NR") else b
            ja_pairs_merged.append((a_m, b_m))

    s2_orig_agree = sum(a == b for a, b in ja_pairs_original) / len(ja_pairs_original)
    s2_orig_k = cohen_kappa(ja_pairs_original)
    s2_orig_ac1 = gwet_ac1(ja_pairs_original)

    s2_merge_agree = sum(a == b for a, b in ja_pairs_merged) / len(ja_pairs_merged)
    s2_merge_k = cohen_kappa(ja_pairs_merged)
    s2_merge_ac1 = gwet_ac1(ja_pairs_merged)

    print("\n=== A5: SENSITIVITY ANALYSIS MERGING N & NR (163 Jointly Applicable Pairs) ===")
    print(f"Original (distinct N vs NR): Agree = {s2_orig_agree*100:.1f}% ({sum(a==b for a,b in ja_pairs_original)}/163), Kappa = {s2_orig_k:.3f}, AC1 = {s2_orig_ac1:.3f}")
    print(f"Merged   (N + NR collapsed): Agree = {s2_merge_agree*100:.1f}% ({sum(a==b for a,b in ja_pairs_merged)}/163), Kappa = {s2_merge_k:.3f}, AC1 = {s2_merge_ac1:.3f}")
    gain = sum(a == b for a, b in ja_pairs_merged) - sum(a == b for a, b in ja_pairs_original)
    print(f"Collapsing N and NR resolves {gain} of the 28 Stage 2 value disagreements!")

    # Output results to JSON
    output_summary = {
        "A2_cluster_bootstrap": {
            "stage_1_point": s1_point,
            "stage_1_ci_lo": s1_ci_lo,
            "stage_1_ci_hi": s1_ci_hi,
            "stage_2_point": s2_point,
            "stage_2_ci_lo": s2_ci_lo,
            "stage_2_ci_hi": s2_ci_hi,
        },
        "A3_stage_1_chance_corrected": {
            "pooled": {
                "raw_agree": pooled_s1_raw,
                "kappa": pooled_s1_k,
                "ac1": pooled_s1_ac1,
            },
            "items": s1_item_results,
        },
        "A4_locator_classification": {
            "total": dict(loc_classes),
            "AD": dict(ad_loc_classes),
            "VD": dict(vd_loc_classes),
        },
        "A5_sensitivity_merged_N_NR": {
            "original_agree_pct": s2_orig_agree * 100,
            "original_kappa": s2_orig_k,
            "original_ac1": s2_orig_ac1,
            "merged_agree_pct": s2_merge_agree * 100,
            "merged_kappa": s2_merge_k,
            "merged_ac1": s2_merge_ac1,
            "disagreements_resolved": gain,
        },
    }

    out_file = REPO_ROOT / "rinah-posthoc-results.json"
    out_file.write_text(json.dumps(output_summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nSaved summary to {out_file}")


if __name__ == "__main__":
    main()
