#!/bin/bash
FASTQ="/home/rothlab/rli/01_ngsdata/210302_MGY360/*_R1*.fastq.gz"
echo $FASTQ
for file in $FASTQ; do
    echo $file
    fname=${file##*/}
    qsub -V -cwd -e mgy_e/$fname".txt" -o mgy_o/$fname".txt" sge_sub.sh $file
done
