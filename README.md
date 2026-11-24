# Persistent lymphopenia and 28-day mortality in ICU patients with pneumonia

R code for a retrospective two-database cohort study (MIMIC-IV v3.1, eICU-CRD v2.0) with a 96-hour landmark design, reported according to STROBE and RECORD.

## Design
- Population: adults, first ICU stay, pneumonia, baseline lymphopenia (first absolute lymphocyte count [ALC] from 6 h before to 24 h after ICU admission <1.0 ×10⁹/L), alive and in hospital at 96 h, at least one ALC during 48–96 h. Excluded: haematological malignancy, HIV/AIDS, solid-organ or stem-cell transplantation.
- Exposure: persistent lymphopenia = last ALC during 48–96 h <1.0 ×10⁹/L (reference: recovery, ≥1.0 ×10⁹/L).
- Primary outcome: all-cause death within 28 days of ICU admission, counted from the 96-h landmark.
- Missing covariates: multiple imputation by chained equations (20 data sets, separately per database); all estimates combined with Rubin's rules.

## Data access
MIMIC-IV and eICU-CRD are available from PhysioNet after credentialing and signing the data use agreement. No patient-level data are included in this repository.

## Run
```
MIMIC_DIR=/path/mimic-iv-3.1 EICU_DIR=/path/eicu-crd-2.0 ./run_all.sh
```
| Script | Purpose |
|---|---|
| R/00_extract.R | Builds a local DuckDB with the MIMIC-IV and eICU-CRD tables used |
| R/00b_extract_lymphocyte.R | Serial lymphocyte/WBC counts with specimen identifiers, diagnoses with sequence numbers, corticosteroids, vasopressors, ventilation, GCS, FiO2, blood cultures |
| R/01_cohort.R | Cohort selection, exposure, covariates (including six-organ first-day SOFA and PaO2/FiO2), outcomes |
| R/02_analysis.R | Multiple imputation; primary, secondary, confounding-control, sensitivity, selection, subgroup and diagnostic analyses; two-stage meta-analysis |
| R/04_figures.R | Figure 1 (flow chart), Figure 2 (Kaplan–Meier), Supplementary Figures S1 (spline) and S2 (subgroups), TIFF 300 dpi |
| R/05_tables.R | Tables 1–3 and Supplementary Tables S2–S6 |

R version 4.6.1. Packages: data.table, DBI, duckdb, mice, survival, splines, WeightIt, MatchIt, metafor, logistf, car, cmprsk, pROC, EValue, sandwich, lmtest, ggplot2, patchwork.

## License
MIT
