#!/usr/bin/env bash
set -euo pipefail

out_dir="reference_genome"
log_dir="results/logs/day1"
assembly="GCF_025998455.1_ASM2599845v1"
base_url="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/025/998/455/${assembly}"

# Fail early if the download host cannot be reached.
echo "Checking network access to NCBI..."
if ! curl -fsS --connect-timeout 10 -o /dev/null https://ftp.ncbi.nlm.nih.gov/; then
  echo "ERROR: NCBI cannot be reached. If curl reports error 6, DNS resolution is failing." >&2
  echo "Current resolver configuration:" >&2
  cat /etc/resolv.conf 2>/dev/null || true
  exit 1
fi

mkdir -p "${out_dir}" "${log_dir}"
log_file="${log_dir}/reference_download.log"
: > "${log_file}"

echo "Downloading genome FASTA" | tee -a "${log_file}"
curl -fL "${base_url}/${assembly}_genomic.fna.gz" -o "${out_dir}/genome.fa.gz"
gzip -df "${out_dir}/genome.fa.gz"

echo "Downloading GFF3 annotation" | tee -a "${log_file}"
curl -fL "${base_url}/${assembly}_genomic.gff.gz" -o "${out_dir}/annotation.gff3.gz"
gzip -df "${out_dir}/annotation.gff3.gz"

echo "Downloading GTF annotation" | tee -a "${log_file}"
curl -fL "${base_url}/${assembly}_genomic.gtf.gz" -o "${out_dir}/annotation.gtf.gz"
gzip -df "${out_dir}/annotation.gtf.gz"

echo "Creating BED12 gene model for RSeQC" | tee -a "${log_file}"
awk -F '\t' '
BEGIN { OFS = "\t" }
$3 == "gene" {
  split($9, attributes, ";")
  sub(/^ID=/, "", attributes[1])
  bed_start = $4 - 1
  print $1, bed_start, $5, attributes[1], 0, $7, bed_start, $5, 0, 1, $5 - bed_start, 0
}
' "${out_dir}/annotation.gff3" > "${out_dir}/genes.bed"

echo "Creating SHA-256 checksums" | tee -a "${log_file}"
(
cd "${out_dir}"
sha256sum genome.fa annotation.gff3 annotation.gtf genes.bed > reference_files.sha256
)

cat > "${out_dir}/REFERENCE.md" <<EOF
# Reference genome and annotation

- Organism: Helicobacter pylori strain CHC155
- RefSeq assembly accession: GCF_025998455.1
- Assembly name: ASM2599845v1
- Chromosome accession: NZ_AP026446.1
- Primary workflow annotation: reference_genome/annotation.gff3
- Companion annotation: reference_genome/annotation.gtf
- RSeQC BED12 gene model: reference_genome/genes.bed
- Genome FASTA: reference_genome/genome.fa
- SHA-256 checksums: reference_genome/reference_files.sha256
- NCBI FTP directory: ${base_url}/
- Download script: scripts/day1/01_reference_genome_annotation.sh
- Downloaded (UTC): $(date -u '+%Y-%m-%dT%H:%M:%SZ')

All workshop reads and all derived coordinates are generated against this
assembly. Processed files from other assemblies are not compatible.
EOF

echo "Reference bundle ready in ${out_dir}" | tee -a "${log_file}"