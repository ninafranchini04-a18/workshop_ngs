#!/usr/bin/env bash
set -euo pipefail

srr="ERR11479311"
sample="WT_2"

outdir="raw_data/fastq/rnaseq"
tmpdir="sra_tmp"

mkdir -p "${outdir}" "${tmpdir}"

fasterq-dump "${srr}" \
  --split-files \
  --threads 2 \
  --progress \
  --details \
  --outdir "${outdir}" \
  --temp "${tmpdir}"

mv "${outdir}/${srr}_1.fastq" "${outdir}/${sample}_R1.fastq"
mv "${outdir}/${srr}_2.fastq" "${outdir}/${sample}_R2.fastq"

gzip "${outdir}/${sample}_R1.fastq" "${outdir}/${sample}_R2.fastq"
