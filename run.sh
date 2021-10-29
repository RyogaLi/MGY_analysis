#!/bin/bash
# path to fastq: /home/rothlab/rli/01_ngsdata/180226_mgy
# output path: ./output/
# bowtie2: /home/rothlab/rli/lib/bowtie2-2.3.4.1-linux-x86_64/bowtie2
file=$1
bowtie2=/home/rothlab/rli/lib/bowtie2-2.3.4.1-linux-x86_64/bowtie2
samtools=/home/rothlab/rli/lib/samtools-1.4.1/samtools
bcftools=/home/rothlab/rli/lib/bcftools-1.4.1/bcftools

FASTQ="/home/rothlab/rli/01_ngsdata/210302_MGY360/*_R1*.fastq.gz"
REF="./ref/S288C_reference_sequence_R64-2-1_20150113"
FASTA="./ref/S288C_reference_sequence_R64-2-1_20150113.fasta"
OUTPUT="./output_21/"
PLOIDY="./sample_ploidy.txt"

#for file in $FASTQ; do
begin=$(date +"%s")

fname=${file##*/} # file basename
fname=${fname%_R1_001.fastq.gz} # filename without extension
#echo "Analyzing.. $fname"
# make a directory for each student
mkdir $OUTPUT$fname

DIR=$OUTPUT$fname
log=$DIR"/"$fname"_log.txt"

# R1 and R2
R1=$file
R2=${file%_R1_001.fastq.gz}"_R2_001.fastq.gz"
# echo $R2
#echo "Aligning...." >> $log
####################################################################
# align with bowtie2
SAM=$DIR"/"$fname"_raw.sam"
# echo $SAM
#bowtie2 -p 20 -x $REF -1 $R1 -2 $R2 -S $SAM 2> $DIR"/align_log.txt"
BAM=$DIR"/"$fname"_raw.bam"
echo "Converting .sam file to .bam file" >> $log
samtools view -S -b $SAM > $BAM
SORTED=$DIR"/"$fname"_sorted"
echo "Sorting .bam file to $SORTED" >> $log
samtools sort $BAM -o $SORTED".bam"
samtools index $SORTED".bam"
###################################################################

echo "Variants call...." >> $log
####################################################################
# variants call
BCF=$DIR"/"$fname"_raw.bcf"
echo "Pile up ...." >> $log
if [[ "$SORTED" =~ 04|05|06 ]]; then
      echo "It's there."
    samtools mpileup -uf $FASTA $SORTED".bam" | bcftools call -vc -> $BCF
else
    samtools mpileup -uf $FASTA $SORTED".bam" | bcftools call --ploidy 1 -vc -> $BCF
fi
VCF=$DIR"/"$fname"_filtered.vcf"
echo "Convert to vcf...." >> $log
bcftools view $BCF | vcfutils.pl varFilter -D1000 > $VCF
    
termin=$(date +"%s")
difftimelps=$(($termin-$begin))
echo "$(($difftimelps / 60)) minutes and $(($difftimelps % 60)) seconds elapsed for Script Execution." >> $log
echo "DONE with $fname" >> $log

####################################################################
    
 #   break

#done

