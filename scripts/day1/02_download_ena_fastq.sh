#!/usr/bin/env bash
set -euo pipefail

mkdir -p \
  raw_data/metadata \
  raw_data/fastq/rnaseq \
  raw_data/fastq/chipseq \
  results/logs/day1

log_file="results/logs/day1/metadata_fastq_download.log"
: > "${log_file}"

# Fail early if EMBL-EBI cannot be reached.
echo "Checking network access to EMBL-EBI..." | tee -a "${log_file}"
if ! curl -fsS --connect-timeout 10 -o /dev/null https://www.ebi.ac.uk/; then
  echo "ERROR: EMBL-EBI cannot be reached. If curl reports error 6, DNS resolution is failing." >&2
  echo "Current resolver configuration:" >&2
  cat /etc/resolv.conf 2>/dev/null || true
  exit 1
fi


# -------------------------------------------------------------------------
# RNA-seq metadata
# -------------------------------------------------------------------------

echo "Downloading RNA-seq metadata" | tee -a "${log_file}"

curl -fL \
  https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/025/E-MTAB-13025/Files/E-MTAB-13025.sdrf.txt \
  -o raw_data/metadata/E-MTAB-13025.sdrf.txt


# Keep only WT samples
grep -v "ΔHP1021" \
  raw_data/metadata/E-MTAB-13025.sdrf.txt \
  > raw_data/metadata/RNAseq_metadata.tmp


# Add two new columns:
# HTTPS_FASTQ_URI = FASTQ address converted from ftp:// to https://
# FASTQ_NAME      = simplified name such as WT_1_R1.fastq.gz

head -n 1 raw_data/metadata/RNAseq_metadata.tmp |
  sed $'s/$/\tHTTPS_FASTQ_URI\tFASTQ_NAME/' \
  > raw_data/metadata/RNAseq_metadata.txt


tail -n +2 raw_data/metadata/RNAseq_metadata.tmp |
while IFS= read -r line; do

  sample=$(printf '%s\n' "${line}" | cut -d $'\t' -f1)
  url=$(printf '%s\n' "${line}" | cut -d $'\t' -f32 | tr -d '\r')

  # ftp://... -> https://...
  url="https://${url#ftp://}"

  # ERR11479310_1.fastq.gz -> 1
  filename="${url##*/}"
  mate="${filename%.fastq.gz}"
  mate="${mate##*_}"

  # WT_1 + 1 -> WT_1_R1.fastq.gz
  fastq_name="${sample}_R${mate}.fastq.gz"

  printf '%s\t%s\t%s\n' \
    "${line}" \
    "${url}" \
    "${fastq_name}" \
    >> raw_data/metadata/RNAseq_metadata.txt

done

rm raw_data/metadata/RNAseq_metadata.tmp


# -------------------------------------------------------------------------
# ChIP-seq metadata
# -------------------------------------------------------------------------

echo "Downloading ChIP-seq metadata" | tee -a "${log_file}"

curl -fL \
  https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/026/E-MTAB-13026/Files/E-MTAB-13026.sdrf.txt \
  -o raw_data/metadata/E-MTAB-13026.sdrf.txt


# Keep only WT samples
grep -v "ΔHP1021" \
  raw_data/metadata/E-MTAB-13026.sdrf.txt \
  > raw_data/metadata/Chipseq_metadata.tmp


# Add the HTTPS_FASTQ_URI and FASTQ_NAME columns

head -n 1 raw_data/metadata/Chipseq_metadata.tmp |
  sed $'s/$/\tHTTPS_FASTQ_URI\tFASTQ_NAME/' \
  > raw_data/metadata/Chipseq_metadata.txt


tail -n +2 raw_data/metadata/Chipseq_metadata.tmp |
while IFS= read -r line; do

  sample=$(printf '%s\n' "${line}" | cut -d $'\t' -f1)
  url=$(printf '%s\n' "${line}" | cut -d $'\t' -f31 | tr -d '\r')

  # ftp://... -> https://...
  url="https://${url#ftp://}"

  # ERR11468761_1.fastq.gz -> 1
  filename="${url##*/}"
  mate="${filename%.fastq.gz}"
  mate="${mate##*_}"

  # WT_1 + 1 -> WT_1_R1.fastq.gz
  fastq_name="${sample}_R${mate}.fastq.gz"

  printf '%s\t%s\t%s\n' \
    "${line}" \
    "${url}" \
    "${fastq_name}" \
    >> raw_data/metadata/Chipseq_metadata.txt

done

rm raw_data/metadata/Chipseq_metadata.tmp


# -------------------------------------------------------------------------
# Download RNA-seq FASTQ files
# -------------------------------------------------------------------------

echo "Downloading RNA-seq FASTQ files" | tee -a "${log_file}"

# The two new columns are columns 33 and 34
rev raw_data/metadata/RNAseq_metadata.txt |
cut -d $'\t' -f1,2 |
rev |
tail -n +2 |
while IFS=$'\t' read -r url fastq_name; do

  echo "Downloading ${fastq_name}" | tee -a "${log_file}"

  curl -fL \
    "${url}" \
    -o "raw_data/fastq/rnaseq/${fastq_name}"

done


# -------------------------------------------------------------------------
# Download ChIP-seq FASTQ files
# -------------------------------------------------------------------------

echo "Downloading ChIP-seq FASTQ files" | tee -a "${log_file}"

# The two new columns are columns 32 and 33
rev raw_data/metadata/Chipseq_metadata.txt |
cut -d $'\t' -f1,2 |
rev |
tail -n +2 |
while IFS=$'\t' read -r url fastq_name; do

  echo "Downloading ${fastq_name}" | tee -a "${log_file}"

  curl -fL \
    "${url}" \
    -o "raw_data/fastq/chipseq/${fastq_name}"

done


# -------------------------------------------------------------------------
# Record what was done
# -------------------------------------------------------------------------

cat > raw_data/metadata/METADATA.md <<EOF
# Metadata download and modification

RNA-seq metadata:
https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/025/E-MTAB-13025/Files/E-MTAB-13025.sdrf.txt

ChIP-seq metadata:
https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/026/E-MTAB-13026/Files/E-MTAB-13026.sdrf.txt

The metadata files were modified as follows:

- Samples carrying the ΔHP1021 genotype were excluded.
- ENA FASTQ addresses were converted from ftp:// to https://.
- A new HTTPS_FASTQ_URI column was added.
- A new FASTQ_NAME column was added.
- FASTQ names combine the sample name with the paired-end read number.

Example:

WT_1 + ERR11479310_1.fastq.gz -> WT_1_R1.fastq.gz
WT_1 + ERR11479310_2.fastq.gz -> WT_1_R2.fastq.gz

RNA-seq FASTQ files:
raw_data/fastq/rnaseq

ChIP-seq FASTQ files:
raw_data/fastq/chipseq
EOF


echo "Metadata processing and FASTQ download complete" | tee -a "${log_file}"