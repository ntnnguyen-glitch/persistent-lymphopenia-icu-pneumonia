#!/usr/bin/env bash
# Usage: MIMIC_DIR=/path/mimic-iv-3.1 EICU_DIR=/path/eicu-crd-2.0 ./run_all.sh
set -euo pipefail
cd "$(dirname "$0")"; mkdir -p output; cd output
export DUCKDB="${DUCKDB:-$PWD/w.duckdb}"
[ -f "$DUCKDB" ] || Rscript ../R/00_extract.R
for s in 00b_extract_lymphocyte 01_cohort 02_analysis 04_figures 05_tables; do Rscript ../R/$s.R > $s.log 2>&1; done
