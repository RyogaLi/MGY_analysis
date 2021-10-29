#!/bin/bash
python=/home/rothlab/rli/py/bin/python2.7
#parental=./output/tao3_parental_1/tao3_parental_1_filtered.vcf
cds=./cds.txt

for file in ./output_21/*/*.vcf; do
  echo $file
  $python compare.py --f $file --c $cds
  echo "OK"
    #break
done
