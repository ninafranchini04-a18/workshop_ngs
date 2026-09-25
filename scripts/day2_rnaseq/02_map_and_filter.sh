#!usr/bin/env bash
set -euo pipefail

path_index="reference_genome/bowtie2_index/genome"
fastq_dir="results/day1_qc/trimmed_fastq/rnaseq"

sam_dir="/Users/ieo5244/workshop_ngs/results/day2_qc/rnaseq/mapped"
log_dir="/Users/ieo5244/workshop_ngs/results/day2_qc/rnaseq/logs"


mkdir -p $sam_dir $log_dir

read1=/Users/ieo5244/workshop_ngs/results/day1_qc/trimmed_fastq/rnaseq/WT_1_R1.trimmed.fastq.gz
read2=/Users/ieo5244/workshop_ngs/results/day1_qc/trimmed_fastq/rnaseq/WT_1_R2.trimmed.fastq.gz

for fastq in ${fastq_dir}/*_R1.*.gz
do 
    echo $fastq
    read1=${fastq}
    read2="${read1/_R1.trimmed.fastq.gz/_R2.trimmed.fastq.gz/}"
    sample="$(basename "{read1}" _R1.trimmed.fastq.gz)"

    bowtie2 \
    -x $path_index \
    -1 $read1 \
    -2 $read2 \
    -p 2 \
    -S ${sam_dir}/${sample}.sam 2> ${log_dir}/${sample}.log

samtools sort \
    @ 2 \
    -o ${sam_dir}/${sample}.sort.bam \
    ${sam_dir}/${sample}.sam



done
  
   
