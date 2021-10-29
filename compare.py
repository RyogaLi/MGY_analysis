import numpy as np
import pandas as pd
import argparse
# compare two vcf file and output the diff

def load(f):
	vcf_file = []
	with open(f, "r") as vcf:
		for line in vcf:
			if not line.startswith("##"):
				vcf_file.append(line.split())
	colnames = vcf_file[0]
	vcf_file = pd.DataFrame(vcf_file[1:], columns=colnames)
	return vcf_file

def compare(f1, f2):
	"""
	compare two files and return a 3rd file with diff
	(use f1 as reference)
	"""
	parental = []
	# parental.append("tao3_parental_1_filtered.vcf")
	for index, row in f2.iterrows():
		# find the variants in f1
		# print row["POS"]
		pa = f1[(f1["#CHROM"]==row["#CHROM"]) & (f1["POS"]==row["POS"]) & (f1["REF"]==row["REF"])]
		if list(pa["ALT"]) == []:
		#	print list(pa["ALT"])
			parental.append("N/A")
		else:
			#print any(i in list(pa["ALT"]) for i in list(row["ALT"]))
			#print list(pa)
			#print list(pa["ALT"])
			if any(i in list(pa["ALT"])[0].split(",") for i in row["ALT"].split(",")):
			    pa = pa["./output/tao3_parental_1/tao3_parental_1_sorted.bam"]
			    parental.append(str(list(pa)[0]))
			else:
		#	    print list(pa["ALT"])[0].split(",")
			    #print row["ALT"].split(",")
			    parental.append("N/A")
			#print row["ALT"]
			#print list(pa["ALT"])
			#parental.append(str(list(pa)[0]))
		# if variant in f1, parental.append(parental col)
	f2["parental"] = parental
	return f2

def find_gene(variants, gene_file):
	genes = pd.read_table(gene_file, sep="\t")
	v_gene = []
	find = 0
	for index, row in variants.iterrows():
		# get genes on that chr
		chromo = genes[genes["Chromosome"]==row["#CHROM"]]
		g = []
		for i, r in chromo.iterrows():
			positions = range(int(r["Start site"]), int(r["End site"]+1))
			v_pos = int(row["POS"])
			if v_pos in positions:
				entry=r["ID"]+"("+r["Orientation"]+")"
				g.append(entry)
				# v_gene.append(r["ID"])
				find = 1
		if find == 0:
			v_gene.append("N/A")
		if find == 1:
			v_gene.append(",".join(g))
		find = 0

	variants["CDS(Orientation)"] = v_gene
	return variants

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Process some integers.')
	parser.add_argument('--p', type=str, help='Parental')
	parser.add_argument('--c', type=str, help='CDS')
	parser.add_argument('--f', type=str, help='Input vcf file')

	args = parser.parse_args()
	#parental = args.p
	vcf = args.f
	cds = args.c

	# load vcf
	#par = load(parental)
	f = load(vcf)

	#final = compare(par, f)
	final = find_gene(f, cds)

	filename = vcf.split("/")[-1].split(".")[0]+"_analyzed.csv"
	final.to_csv("./csv_files2021/"+filename)
