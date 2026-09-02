
# Interrater calibration materials

This repository supports the manuscript "Applicability, not judgement, drives disagreement when auditing safety reporting in nursing generative AI studies: a methodological calibration study."

Two independent human reviewers coded the reports. The public files label their paired outputs as coder A and coder B. Adjudicated values were confirmed by two reviewers.

## Contents

- `data/coding-dictionary.csv`: candidate 70-item coding dictionary.
- `data/paired-coding-values.csv`: 312 paired values from the independent 12-paper by 26-item calibration sample, labeled only as coder A and coder B.
- `data/adjudicated-disagreements.csv`: final values for the 75 disagreements, with confirmation by two reviewers.
- `data/rc5-reliability.csv`: item-level agreement, coefficients, Wilson intervals, and bootstrap intervals.
- `data/rc6-item-disposition.csv`: item-level post-calibration analytic classifications.
- `data/rc6-tier-summary.csv`: counts across the five analytic classifications.
- `scripts/`: dependency-free Python code used for agreement coefficients and uncertainty intervals.

## Reproduction

From the repository root, run:

```bash
python3 scripts/calibration_intervals.py
```

The command reports the number of analyzed items and the result of the stricter dual-coefficient sensitivity analysis. Python 3.10 or later is recommended; no third-party Python packages are required for this calculation.

## Release boundary

The paired dataset excludes evidence locators, free-text reviewer notes, and full-text excerpts. Copyrighted complete texts and restricted source documents are not distributed. Paper identifiers are retained because the unit of analysis is the published report and source identity is required for auditability. All released coding values underwent final review by two researchers.

## License

- Data and documentation: CC BY 4.0 (`LICENSE-DATA.md`).
- Source code: MIT (`LICENSE-CODE.md`).

## Release status

Version 0.1.2 is available at https://github.com/871922384/nursing-genai-audit-calibration. The versioned release is https://github.com/871922384/nursing-genai-audit-calibration/releases/tag/v0.1.2. Paired coding values are unchanged from v0.1.0. This release aligns post-calibration analytic classifications with the manuscript primary rule after the jointly-applicable and bootstrap display filters: 2 of 26 items met the primary rule, with analytic-tier counts 44/2/10/12/2. Documentation names the current manuscript title, two independent human reviewers, and the retrospective OSF Open-Ended Registration DOI 10.17605/OSF.IO/5KQD9 (https://osf.io/5kqd9/; associated project https://osf.io/kes4f/; registered 2 September 2026). The OSF registration postdates data collection.
