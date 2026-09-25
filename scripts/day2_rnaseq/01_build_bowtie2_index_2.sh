#!/usr/bin/env bash
set -euo pipefail

genome="reference_genome/genome.fa"
index_path="reference_genome/bowtie2_index/genome"

mkdir -p reference_genome/bowtie2_index

bowtie2-build --threads 2 ${genome} ${index_path}