#!/usr/bin/env bash
set -euo pipefail

genome="reference_genome/genome.fa"
index_prefix="reference_genome/bowtie2_index/genome"

mkdir -p reference_genome/bowtie2_index

if [[ ! -f "${genome}" ]]; then
  echo "ERROR: genome FASTA not found: ${genome}" >&2
  exit 1
fi

bowtie2-build \
  --threads 2 \
  "${genome}" \
  "${index_prefix}"