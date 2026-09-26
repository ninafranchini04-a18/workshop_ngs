#!/usr/bin/env bash
set -euo pipefail

bam_dir="results/day2_rnaseq/bam"
out_dir="results/day2_rnaseq/bigwig"
scale_file="results/day2_rnaseq/deseq2/deseq2_size_factors.tsv"

mkdir -p "${out_dir}"

if [[ ! -f "${scale_file}" ]]; then
  echo "ERROR: DESeq2 scale-factor table not found: ${scale_file}" >&2
  exit 1
fi

tail -n +2 "${scale_file}" |
while IFS=$'\t' read -r sample deseq2_size_factor bamcoverage_scale_factor
do
  bam="${bam_dir}/${sample}.filtered.bam"

  if [[ ! -f "${bam}" ]]; then
    echo "ERROR: BAM file not found: ${bam}" >&2
    exit 1
  fi

  echo "${sample}: DESeq2 size factor=${deseq2_size_factor}; bamCoverage scale factor=${bamcoverage_scale_factor}"

bamCoverage \
    -b "${bam}" \
    --filterRNAstrand forward \
    -o "${out_dir}/${sample}.norm.fw.bw" \
    --scaleFactor "${bamcoverage_scale_factor}" \
    --binSize 10 \
    --numberOfProcessors 2
 bamCoverage \
    -b "${bam}" \
    --filterRNAstrand reverse \
    -o "${out_dir}/${sample}.norm.rv.bw" \
    --scaleFactor "${bamcoverage_scale_factor}" \
    --binSize 10 \
    --numberOfProcessors 2
done