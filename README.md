# Reproducible NGS workshop
## Day 1

Day 1 created the reproducible project structure, downloaded the matched reference genome and annotation, converted the RNA-seq and ChIP-seq SDRFs into the two metadata tables used by the workshop, downloaded and renamed the WT FASTQ files, and used fastp to preprocess and assess the reads.

After each sample was processed successfully, its raw FASTQ pair was removed to save space. Large files such as trimmed FASTQ files and reference sequence/annotation files are ignored by Git and are not pushed to GitHub. Small provenance records and quality reports are committed.