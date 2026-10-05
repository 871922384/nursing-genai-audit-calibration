
# Interrater calibration materials

This repository supports the manuscript "Why Reviewers Disagree When Auditing Reporting in Nursing Generative AI Studies: An Applicability-First Interrater Agreement Study"

Two independent human reviewers coded the reports. The public files label their paired outputs as coder A and coder B. Adjudicated values were confirmed by two reviewers.

## Contents

- `data/coding-dictionary.csv`: candidate 70-item coding dictionary.
- `data/item-id-map.csv`: canonical item mapping (A01–A26) for the 26 calibration variables.
- `data/paired-coding-values.csv`: 312 paired values from the independent 12-paper by 26-item calibration sample, labeled only as coder A and coder B.
- `data/bootstrap-paper-order.csv`: the 12-paper input sequence used for the published item-level bootstrap intervals. The script applies this sequence even if the paired-value CSV is reordered.
- `data/pilot-screening-exclusions.csv`: 25 distinct preliminary feasibility reports used before the formal 20-report development phase; these reports were excluded from both the development counts and the independent 12-report calibration sample.
- `data/adjudicated-disagreements.csv`: final values for the 75 disagreements, with confirmation by two reviewers.
- `data/disagreement-locators.csv`: sanitized coder A and coder B evidence-locator strings and discrepancy classifications for all 75 disagreements.
- `data/rc5-reliability.csv`: item-level agreement, coefficients, Wilson intervals, and bootstrap intervals.
- `data/rc6-item-disposition.csv`: item-level post-calibration analytic classifications.
- `data/rc6-tier-summary.csv`: counts across the five analytic classifications.
- `scripts/calibration_intervals.py`: dependency-free Python code used for agreement coefficients and uncertainty intervals.
- `scripts/reliability_metrics.py`: recomputes the 26 item-level status and numerical reliability fields from the released paired values.
- `scripts/compute_posthoc_analyses.py`: Python code for the four post hoc analyses (cluster bootstrap, pooled Stage 1 AC1, locator discrepancy classification, N/NR collapsing).

## Reproduction

From the repository root, run:

```bash
python3 scripts/calibration_intervals.py
python3 scripts/reliability_metrics.py --paired data/paired-coding-values.csv --dictionary data/coding-dictionary.csv --item-map data/item-id-map.csv --output reproduced-rc5.csv --verify-against data/rc5-reliability.csv
python3 scripts/compute_posthoc_analyses.py
```

Python 3.10 or later is recommended; no third-party Python packages are required for these calculations.

The item-level bootstrap uses 10,000 resamples, a variable-specific SHA-256 seed, and the public `bootstrap-paper-order.csv` sequence. `calibration_intervals.py` now checks all 234 computed interval and denominator fields against the released RC5 table. The reliability command checks all 17 base fields for each of the 26 items. Five sparse items have no jointly applicable report pair; the public status rule assigns `insufficient_applicable_sample` before applying the agreement thresholds. Coefficients calculated from applicability-discordant pairs remain visible, but do not override this status rule.

## Release boundary

The release includes sanitized reviewer evidence-locator strings and discrepancy classifications. Full-text excerpts, free-text coding notes, complete copyrighted texts, and restricted source documents are not distributed. Paper identifiers are retained because the unit of analysis is the published report and source identity is required for auditability. The published paired coder A and coder B values are frozen pre-adjudication ratings; consensus adjudication did not overwrite either reviewer's independent values. The separate adjudicated-disagreements file records the resolved values.

## License

- Data and documentation: CC BY 4.0 (`LICENSE-DATA.md`).
- Source code: MIT (`LICENSE-CODE.md`).

## Archival citation

The matching version 1.0.3 data and code archive is deposited at Zenodo: https://doi.org/10.5281/zenodo.23156507. Data and documentation use CC BY 4.0; Python source code uses MIT. The OSF DOI below identifies the retrospective registration separately.

## Release status

Submission-linked release v1.0.3 is available at https://github.com/871922384/nursing-genai-audit-calibration/releases/tag/v1.0.3. This reproducibility correction publishes the original bootstrap input order and encodes the sparse-item status rule so that the released paired data reproduce the published RC5 fields. The 312 paired values, statistical results, 75 adjudicated disagreements, and final analytic classifications are unchanged from v1.0.2. Documentation names the current manuscript title, two independent human reviewers, and the retrospective OSF Open-Ended Registration DOI 10.17605/OSF.IO/5KQD9 (https://osf.io/5kqd9/; associated project https://osf.io/kes4f/; registered 2 September 2026). The OSF registration postdates data collection. Development-round counts in the manuscript are descriptive formative records; development-phase paired records are not included in this repository.
