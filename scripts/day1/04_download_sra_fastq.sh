#!/usr/bin/env bash
set -euo pipefail

root="raw_data_sra"

mkdir -p \
  "${root}/metadata" \
  "${root}/sra" \
  "${root}/fastq/rnaseq" \
  "${root}/fastq/chipseq"

# Metadata come from EMBL-EBI; sequencing runs are retrieved through SRA.
echo "Checking network access to EMBL-EBI and NCBI..."
for url in https://www.ebi.ac.uk/ https://www.ncbi.nlm.nih.gov/; do
  if ! curl -fsS --connect-timeout 10 -o /dev/null "${url}"; then
    echo "ERROR: ${url} cannot be reached. If curl reports error 6, DNS resolution is failing." >&2
    echo "Current resolver configuration:" >&2
    cat /etc/resolv.conf 2>/dev/null || true
    exit 1
  fi
done


# -------------------------------------------------------------------------
# RNA-seq metadata
# -------------------------------------------------------------------------

rnaseq_study="ERP147707"

curl -fL \
  "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${rnaseq_study}&result=read_run&fields=run_accession,sample_title,library_layout&format=tsv" \
  -o "${root}/metadata/RNAseq_metadata_all.tsv"

# Keep WT_1-3 and WTS_1-3 and add the local FASTQ prefix
awk -F $'\t' 'BEGIN {OFS="\t"}
NR==1 {print $1,$2,$3,"FASTQ_PREFIX"; next}
$2 ~ /^(WT|WTS)_[123]$/ {print $1,$2,$3,$2}
' "${root}/metadata/RNAseq_metadata_all.tsv" \
  > "${root}/metadata/RNAseq_metadata.txt"


# -------------------------------------------------------------------------
# ChIP-seq metadata
# -------------------------------------------------------------------------

chipseq_study="ERP147583"

curl -fL \
  "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${chipseq_study}&result=read_run&fields=run_accession,sample_title,library_layout&format=tsv" \
  -o "${root}/metadata/Chipseq_metadata_all.tsv"

# Keep WT_1-3 and WTS_1-3 and add the local FASTQ prefix
awk -F $'\t' 'BEGIN {OFS="\t"}
NR==1 {print $1,$2,$3,"FASTQ_PREFIX"; next}
$2 ~ /^(WT|WTS)_[123]$/ {print $1,$2,$3,$2}
' "${root}/metadata/Chipseq_metadata_all.tsv" \
  > "${root}/metadata/Chipseq_metadata.txt"


# -------------------------------------------------------------------------
# Download RNA-seq FASTQ files
# -------------------------------------------------------------------------

tail -n +2 "${root}/metadata/RNAseq_metadata.txt" |
while IFS=$'\t' read -r run sample layout prefix; do

  echo "Downloading RNA-seq ${sample} (${run}); layout=${layout}"

  # Download the SRA run first. If the transfer is interrupted,
  # rerunning prefetch can resume the existing download.
  prefetch --max-size u -O "${root}/sra" "${run}"

  # Validate the local archive object before converting it to FASTQ.
  vdb-validate "${root}/sra/${run}"

  if [[ "${layout}" == "PAIRED" ]]; then
    fastq-dump \
      --split-files \
      --skip-technical \
      --outdir "${root}/fastq/rnaseq" \
      "${root}/sra/${run}"

    mv "${root}/fastq/rnaseq/${run}_1.fastq" "${root}/fastq/rnaseq/${prefix}_R1.fastq"
    mv "${root}/fastq/rnaseq/${run}_2.fastq" "${root}/fastq/rnaseq/${prefix}_R2.fastq"

    gzip "${root}/fastq/rnaseq/${prefix}_R1.fastq"
    gzip "${root}/fastq/rnaseq/${prefix}_R2.fastq"

  else
    fastq-dump \
      --skip-technical \
      --outdir "${root}/fastq/rnaseq" \
      "${root}/sra/${run}"

    mv "${root}/fastq/rnaseq/${run}.fastq" "${root}/fastq/rnaseq/${prefix}_R1.fastq"
    gzip "${root}/fastq/rnaseq/${prefix}_R1.fastq"
  fi

done


# -------------------------------------------------------------------------
# Download ChIP-seq FASTQ files
# -------------------------------------------------------------------------

tail -n +2 "${root}/metadata/Chipseq_metadata.txt" |
while IFS=$'\t' read -r run sample layout prefix; do

  echo "Downloading ChIP-seq ${sample} (${run}); layout=${layout}"

  # Download the SRA run first. If the transfer is interrupted,
  # rerunning prefetch can resume the existing download.
  prefetch --max-size u -O "${root}/sra" "${run}"

  # Validate the local archive object before converting it to FASTQ.
  vdb-validate "${root}/sra/${run}"

  if [[ "${layout}" == "PAIRED" ]]; then
    fastq-dump \
      --split-files \
      --skip-technical \
      --outdir "${root}/fastq/chipseq" \
      "${root}/sra/${run}"

    mv "${root}/fastq/chipseq/${run}_1.fastq" "${root}/fastq/chipseq/${prefix}_R1.fastq"
    mv "${root}/fastq/chipseq/${run}_2.fastq" "${root}/fastq/chipseq/${prefix}_R2.fastq"

    gzip "${root}/fastq/chipseq/${prefix}_R1.fastq"
    gzip "${root}/fastq/chipseq/${prefix}_R2.fastq"

  else
    fastq-dump \
      --skip-technical \
      --outdir "${root}/fastq/chipseq" \
      "${root}/sra/${run}"

    mv "${root}/fastq/chipseq/${run}.fastq" "${root}/fastq/chipseq/${prefix}_R1.fastq"
    gzip "${root}/fastq/chipseq/${prefix}_R1.fastq"
  fi

done
