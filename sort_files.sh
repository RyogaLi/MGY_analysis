#FASTQ="/home/rothlab/rli/www/html/MGY_2019/190219_fastq_20k/"
#CSV="/home/rothlab/rli/02_dev/mgy_pipeline/csv_files/"
OUTPUT="/home/rothlab/rli/02_dev/mgy_pipeline/output_21/"

for f in $(find "$1" -mindepth 1 -maxdepth 1 -type d); do
  echo $f
    zip -r $f.zip $f
    #NAME=$(basename $f)
  #ID=$(echo $NAME | cut -d"_" -f1)
  #IDS=$(echo $NAME | cut -d"_" -f2)
  #BASENAME=$ID\_$IDS
  #DIR=$OUTPUT$BASENAME/
  # move csv files
  #cp $f $DIR
  # move fastq files
  #R1=$FASTQ$BASENAME\_R1.fastq
  #R2=$FASTQ$BASENAME\_R2.fastq
  #cp $R1 $DIR
  #cp $R2 $DIR
  #rm $DIR*.fastq
  #cd $DIR
  #rm *_raw.bam
  #rm *_raw.sam
  #cd $DIR
  #rm *_log.txt
  #cd $(dirname $DIR)
  #echo $BASENAME
  #zip -r $BASENAME.zip $BASENAME/
done
