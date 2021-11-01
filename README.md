#### How to run MGY analysis

0. Pre-built reference file located on galen: `/home/rothlab/rli/02_dev/mgy_pipeline/ref/S288C_reference_sequence_R64-2-1_20150113`

```
usage: main.py [-h] [--fastq FASTQ] [--output OUTPUT] [--ref REF]
               [--ploidy PLOIDY]

MGY analysis

optional arguments:
  -h, --help       show this help message and exit
  --fastq FASTQ    Path to all fastq files you want to analyze
  --output OUTPUT  Output path for sam files
  --ref REF        path to all reference files
  --ploidy PLOIDY  Ploidity (default set to 1)
```

1. Run the alignment and variant call (submit jobs to slurm). Command example:
`python main.py --fastq path/to/fastq --output path/to/output --ref path/to/ref`

2. Run the analysis script for individual vcf file (no jobs submitted, takes around 1-2hrs)
```
usage: compare.py [-h] [--p P] [--c C] [--f F]

Process some integers.

optional arguments:
  -h, --help  show this help message and exit
  --c C       CDS
  --f F       Input vcf file
```
Example
`python compare.py --c path/to/cds.txt --f output/from/main`
