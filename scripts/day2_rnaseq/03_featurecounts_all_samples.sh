#!/usr/bin/env bash
set -euo pipefail

annotation="reference_genome/annotation.gff3"
bam_dir="results/day2_rnaseq/bam"
counts_dir="results/day2_rnaseq/counts"
log_dir="results/logs/day2"
featurecounts_out="${counts_dir}/featureCounts_all_samples.txt"

mkdir -p "${counts_dir}" "${log_dir}"

if [[ ! -f "${annotation}" ]]; then
  echo "ERROR: annotation file not found: ${annotation}" >&2
  exit 1
fi

bam_files=("${bam_dir}"/*.filtered.bam)

if [[ ! -e "${bam_files[0]}" ]]; then
  echo "ERROR: No filtered BAM files found in ${bam_dir}" >&2
  exit 1
fi

echo "Counting reads across ${#bam_files[@]} samples"

featureCounts \
  -T 2 \
  -p \
  --countReadPairs \
  -B \
  -C \
  -s 2 \
  -t gene \
  -g ID \
  -a "${annotation}" \
  -o "${featurecounts_out}" \
  "${bam_files[@]}" \
  2>&1 | tee "${log_dir}/featureCounts_all_samples.log"

echo "Wrote ${featurecounts_out}"
