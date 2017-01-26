#!/bin/bash

esearch -db sra -query $1 | efetch --format runinfo | cut -d ',' -f 1 | grep SRR >srr.txt;
cat srr.txt | while read LINE
do
SUB=$(echo $LINE | cut -c 1-6)
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/$SUB/$LINE/$LINE.sra
fastq-dump --skip-technical --split-files --readids --dumpbase --clip --ncbi_error_report never --log-level fatal $LINE.sra
mv $LINE*1.fastq "`basename $LINE*1.fastq 1.fastq`forward.fq"
mv $LINE*2.fastq "`basename $LINE*2.fastq 2.fastq`reverse.fq"
done
rm *.sra


