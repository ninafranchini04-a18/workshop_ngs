# Workflow execution log

## Day 1

- Created `.gitignore` and `README.md`.
- Connected the local repository to GitHub.
- Tested Docker with a mounted project folder.
- Downloaded the matched `GCF_025998455.1` FASTA, GFF3, and GTF files for chromosome `NZ_AP026446.1` and derived the RSeQC BED12 gene model from the GFF3.
- Ran `02_download_ena_fastq.sh` to download and filter the original SDRFs, create `RNAseq_metadata.txt` and `Chipseq_metadata.txt`, and retrieve the WT FASTQ files.
- Recorded the metadata sources and transformations in `raw_data/metadata/METADATA.md` and reviewed GEO/SRA as an optional route.
- Ran fastp to trim adapters and low-quality tails, filter poor reads, and create trimmed FASTQ files.
- Inspected the fastp HTML reports; raw pairs were removed automatically after successful processing.