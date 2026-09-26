#!/usr/bin/env bash
set -euo pipefail

fastq_dir="results/day1_qc/trimmed_fastq/rnaseq"
index_prefix="reference_genome/bowtie2_index/genome"
bam_dir="results/day2_rnaseq/bam"
log_dir="results/logs/day2"

mkdir -p "${bam_dir}" "${log_dir}"

for read1 in "${fastq_dir}"/*_R1.trimmed.fastq.gz
do
  [[ -e "${read1}" ]] || continue

  read2="${read1/_R1.trimmed.fastq.gz/_R2.trimmed.fastq.gz}"
  sample_id="$(basename "${read1}" _R1.trimmed.fastq.gz)"
  sam_file="${bam_dir}/${sample_id}.sam"
  bam_file="${bam_dir}/${sample_id}.sorted.bam"
  filtered_bam="${bam_dir}/${sample_id}.filtered.bam"

  if [[ ! -f "${read2}" ]]; then
    echo "ERROR: paired read not found for ${read1}" >&2
    exit 1
  fi

  echo "Mapping ${sample_id}"
  bowtie2 \
    -x "${index_prefix}" \
    -1 "${read1}" \
    -2 "${read2}" \
    -p 2 \
    -S "${sam_file}" \
    2> "${log_dir}/${sample_id}.bowtie2.log"

  samtools sort \
    -@ 2 \
    -o "${bam_file}" \
    "${sam_file}"

  samtools index "${bam_file}"

  samtools flagstat "${bam_file}" > "${log_dir}/${sample_id}.flagstat.txt"

  samtools view \
    -@ 2 \
    -b \
    -q 10 \
    -f 2 \
    -F 2828 \
    -o "${filtered_bam}" \
    "${bam_file}"

  samtools index "${filtered_bam}"

  samtools flagstat "${filtered_bam}" > "${log_dir}/${sample_id}.filtered.flagstat.txt"

  rm -f "${sam_file}"
done