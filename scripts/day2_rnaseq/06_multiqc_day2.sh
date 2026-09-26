#!/usr/bin/env bash
set -euo pipefail

day1_fastp_dir="results/logs/day1/fastp"
day2_log_dir="results/logs/day2"
counts_dir="results/day2_rnaseq/counts"
qc_dir="results/day2_rnaseq/qc"

mkdir -p "${qc_dir}"

multiqc \
  "${day1_fastp_dir}" \
  "${day2_log_dir}" \
  "${counts_dir}" \
  --outdir "${qc_dir}" \
  --filename multiqc_report.html \
  --title "Workshop RNA-seq QC: Day 1 + Day 2" \
  --force