#!/usr/bin/env bash
set -euo pipefail

raw_dir="raw_data/fastq"
out_dir="results/day1_qc/trimmed_fastq"
report_dir="results/logs/day1/fastp"

mkdir -p "${out_dir}/rnaseq" "${out_dir}/chipseq" "${report_dir}"

for assay in rnaseq chipseq
do
  for read1 in "${raw_dir}/${assay}"/*_R1.fastq.gz
  do
    [[ -e "${read1}" ]] || continue

    read2="${read1/_R1.fastq.gz/_R2.fastq.gz}"
    sample_id="$(basename "${read1}" _R1.fastq.gz)"

    if [[ ! -f "${read2}" ]]; then
      echo "ERROR: paired read not found for ${read1}" >&2
      exit 1
    fi

    echo "Trimming ${assay} ${sample_id}"
    fastp \
      -i "${read1}" \
      -I "${read2}" \
      -o "${out_dir}/${assay}/${sample_id}_R1.trimmed.fastq.gz" \
      -O "${out_dir}/${assay}/${sample_id}_R2.trimmed.fastq.gz" \
      --detect_adapter_for_pe \
      --cut_tail \
      --cut_tail_window_size 4 \
      --cut_tail_mean_quality 20 \
      --thread 2 \
      --html "${report_dir}/${assay}_${sample_id}.fastp.html" \
      --json "${report_dir}/${assay}_${sample_id}.fastp.json" \
      > "${report_dir}/${assay}_${sample_id}.log" 2>&1

    rm "${read1}" "${read2}"
    echo "Finished ${assay} ${sample_id}; raw FASTQ pair removed"
  done
done