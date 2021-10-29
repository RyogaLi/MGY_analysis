#!/usr/bin/env python3.7
import argparse
import os
import glob


def main(arguments):
    """
    Submit alignment and variant call jobs
    """
    # make output dir if it doesn't exist
    if not os.path.isdir(arguments.output):
        os.mkdir(arguments.output)
    # go through all the fastq files
    # each fastq file is for one student sample
    file_list = glob.glob(f"{arguments.fastq}/*.fastq.gz")
    fasta = arguments.ref + ".fasta"
    for f in file_list:
        # get Read 1
        if "_R1" not in f: continue
        # in the output dir, make sh dir to store all the sh files
        sh_dir = os.path.join(arguments.output, "GALEN_jobs")
        if not os.path.isdir(sh_dir):
            os.mkdir(sh_dir)

        r2 = f.replace("_R1", "_R2")
        basename = os.path.basename(f)
        sample_id = os.path.basename(f).split("_")[0]
        # make output for this sample
        output_dir = os.path.join(arguments.output, f"{sample_id}")
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)
        sh_file = os.path.join(sh_dir, f"{basename}.sh")

        sam_file = os.path.join(output_dir, basename.replace(".fastq.gz", ".sam"))

        # alignment command
        cmd = f"bowtie2 -p 8 -x {arguments.ref} -1 {f} -2 {r2} -S {sam_file} 2> {output_dir}/align_log.txt\n"
        # convert to bam file
        bam_file = sam_file.replace(".sam", ".bam")
        cmd += f"samtools view -S -b {sam_file} > {bam_file}\n"
        # sort and index bam file
        bam_sorted = bam_file.replace(".bam", "_sorted.bam")
        cmd += f"samtools sort {bam_file} -o {bam_sorted}\n"
        cmd += f"samtools index {bam_sorted} {bam_sorted.replace('.bam', '.bai')}\n"
        # variant call
        bcf = sam_file.replace(".sam", ".bcf")
        cmd += f"samtools mpileup -f {fasta} {bam_sorted} | bcftools call --ploidy {arguments.ploidy} -vc -> {bcf}\n"
        # convert to vcf
        vcf = bcf.replace(".bcf", ".vcf")
        #cmd += f"bcftools view {bcf} | vcfutils.pl varFilter -D1000 > {vcf}"
        #cmd += f"bcftools mpileup --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,"
        #             f"INFO/ADF,INFO/ADR -f {fasta} {bam_sorted)} > {bcf}\n"
        cmd += f"bcftools view -u {bcf} > {vcf}\n"
        # write header to sh file
        header = f"#!/bin/bash\n#SBATCH --time=24:00:00\n#SBATCH --job-name={basename}\n#SBATCH " \
             f"--cpus-per-task=8\n#SBATCH --error={sam_file.replace('.sam', '')}-%j.log\n#SBATCH --mem=10G\n#SBATCH " \
             f"--output={sam_file.replace('.sam', '')}-%j.log\n"

        with open(sh_file, "w") as sh_f:
            sh_f.write(header)
            sh_f.write(cmd)
        os.system(f"sbatch {sh_file}")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='MGY analysis')

    # parameters for cluster
    parser.add_argument("--fastq", help="Path to all fastq files you want to analyze")
    parser.add_argument("--output", help="Output path for sam files")
    parser.add_argument("--ref", help="path to all reference files")
    parser.add_argument("--ploidy", type=int, help="Ploidity", default=1)
    args = parser.parse_args()
    main(args)
