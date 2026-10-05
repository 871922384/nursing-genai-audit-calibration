
# Interrater calibration materials

This repository supports the manuscript "Why Reviewers Disagree When Auditing Reporting in Nursing Generative AI Studies: An Applicability-First Interrater Agreement Study"

Two independent human reviewers coded the reports. The public files label their paired outputs as coder A and coder B. Adjudicated values were confirmed by two reviewers.

## Contents

- `data/coding-dictionary.csv`: candidate 70-item coding dictionary.
- `data/item-id-map.csv`: canonical item mapping (A01–A26) for the 26 calibration variables.
- `data/paired-coding-values.csv`: 312 paired values from the independent 12-paper by 26-item calibration sample, labeled only as coder A and coder B.
- `data/adjudicated-disagreements.csv`: final values for the 75 disagreements, with confirmation by two reviewers.
- `data/disagreement-locators.csv`: sanitized coder A and coder B evidence-locator strings and discrepancy classifications for all 75 disagreements.
- `data/rc5-reliability.csv`: item-level agreement, coefficients, Wilson intervals, and bootstrap intervals.
- `data/rc6-item-disposition.csv`: item-level post-calibration analytic classifications.
- `data/rc6-tier-summary.csv`: counts across the five analytic classifications.
- `scripts/calibration_intervals.py`: dependency-free Python code used for agreement coefficients and uncertainty intervals.
- `scripts/compute_posthoc_analyses.py`: Python code for the four post hoc analyses (cluster bootstrap, pooled Stage 1 AC1, locator discrepancy classification, N/NR collapsing).

## Reproduction

From the repository root, run:

```bash
python3 scripts/calibration_intervals.py
python3 scripts/compute_posthoc_analyses.py
```

Python 3.10 or later is recommended; no third-party Python packages are required for these calculations.

## Release boundary

The release includes sanitized reviewer evidence-locator strings and discrepancy classifications. Full-text excerpts, free-text coding notes, complete copyrighted texts, and restricted source documents are not distributed. Paper identifiers are retained because the unit of analysis is the published report and source identity is required for auditability. The published paired coder A and coder B values are frozen pre-adjudication ratings; consensus adjudication did not overwrite either reviewer's independent values. The separate adjudicated-disagreements file records the resolved values.

## License

- Data and documentation: CC BY 4.0 (`LICENSE-DATA.md`).
- Source code: MIT (`LICENSE-CODE.md`).

## Release status

Submission-linked release v1.0.2 is available at https://github.com/871922384/nursing-genai-audit-calibration/releases/tag/v1.0.2. This documentation-only correction fixes the package-version identifier in RELEASE-STATUS.md and clarifies the status of independent paired values. The release includes the canonical item identifier map (A01–A26), sanitized reviewer locator strings, discrepancy classifications, and post hoc computation scripts. The independent calibration data and analytic results are unchanged from v1.0.1. Documentation names the current manuscript title, two independent human reviewers, and the retrospective OSF Open-Ended Registration DOI 10.17605/OSF.IO/5KQD9 (https://osf.io/5kqd9/; associated project https://osf.io/kes4f/; registered 2 September 2026). The OSF registration postdates data collection. Development-round counts in the manuscript are descriptive formative records; development-phase paired records are not included in this repository.
