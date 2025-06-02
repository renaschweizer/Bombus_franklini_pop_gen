#!/usr/bin/env python3

'''
This script creates a consensus diploid sequence in the form of a fastq file. Two
files are read in: a reference sequence (fasta) and a .vcf file for a single individual.
The vcf file is scanned for sites that are heterozygous within the individual,
and then inserts the appropriate IUPAC ambiguity code into that position of the 
reference sequence; all invariant sites are simply the reference base/allele. 
The base quality scores associated with the sequence are all arbitrarily given 
the best score ("~"). This fastq is then used by the f2psmcfa utility to generate
the input file for PSMC analysis. A .csv file with individual heterozygosity
measurements is also written out.
'''

#Example usage
#python make_consensus_diploid_sequence_from_vcf.py msprime_model_origin_100k_shrink_50k_pop_5k_one_ind.vcf msprime_model_origin_100k_shrink_50k_pop_5k.fq msprime_model_origin_100k_shrink_50k_pop_5k.csv

from Bio import SeqIO
import os
import pandas as pd
from pysam import VariantFile
import sys

# Path to the reference genome
reference_file = '/home/rena.schweizer/scripts/simulation_pipeline/B_affinis_chr_1_13mbp.fna'
ref_dict = SeqIO.to_dict(SeqIO.parse(reference_file, "fasta"))

# Number and length of chromosomes
num_chrs = 1
chrom_length = int(1.3e7)

# Paths to input .vcf file, output fastq file, and .csv with individual heterzygosity estimates
vcf_file = sys.argv[1]
output_fq_file = sys.argv[2]
output_het_file = sys.argv[3]

#A dictionary where keys are the two bases and values are the corresponding IUPAC ambiguity code
#"01" is for the Hu et al. vcf output with 0s and 1s rather than bases
amb_codes = {"01": "Y", "AG": "R", "GT": "K", "CG": "S", "CT": "Y", "AC": "M", "AT": "W"}

def vcf_to_psmcfa(vcf_file, output_fq_file, output_het_file):
	
	#Read in the vcf file
	vcf_in = VariantFile(filename=vcf_file)
	
	#Get the individual ID
	sample_name = vcf_in.header.samples[0]
	
	#Make a string of fake quality scores
	qual_scores_list = ['~'] * chrom_length
	qual_scores = ''.join(qual_scores_list)
	
	#Count the number of variable sites to calculate per-individual heterozygosity
	num_var_sites = 0
	
	#Do the conversion one chromosome at a time
	for chrom_num in range(1, (num_chrs + 1)):

		for entry in ref_dict.values():
			#print("Now processing chromosome", chrom_num)	
			ref_seq = entry.seq
			#Find the variant positions
			for record in vcf_in.fetch():
				chrom = record.contig
				ref_allele = record.ref
				alt_allele = record.alts[0]
				#Get the positions for the specific chromosome
				if str(chrom_num) == str(chrom):				
					position = record.pos - 1 # VCF is 1-based, FASTA is 0-based
					num_alt_alleles = record.info["AC"][0] #get the count of alternate alleles
					if num_alt_alleles == 1: #individual is heterozygous, 0 = homozygous for reference allele, 2 = homozygous for alternate allele
						num_var_sites = num_var_sites + 1
						#Find which IUPAC ambiguity the heterozygote is based on the reference and alternate alleles
						both_bases = ''.join(sorted(ref_allele + alt_allele))
						if both_bases in amb_codes:
							amb = amb_codes[both_bases] #Get the ambiguity code
						ref_seq = ref_seq[:position] + str(amb) + ref_seq[position + 1:] #Put the ambiguity code in the reference sequence
			
			#Write out the new sequence and quality scores to the fastq output
			with open(output_fq_file, 'a') as out:
				out.write(f"@{chrom_num}\n")
				out.write(f"{ref_seq}\n+\n")
				out.write(f"{qual_scores}\n")
				
	#Calculate the individual's total heterozygosity
	gen_wide_het = round((num_var_sites / (num_chrs * chrom_length)), 5)
	
	#Make a pandas data frame with individual ID and heterozygosity
	het_df = pd.DataFrame(columns=['Individual', 'Heterozygosity'])
	het_df.loc[1] = [str(sample_name), gen_wide_het]
	
	#Write the heterozygosity information to file
	if os.path.exists(output_het_file): #het. file exists, append to it
		het_df.to_csv(output_het_file, index = False, header = False, mode = "a")
	
	else: #het. file doesn't exist, write information with header
		het_df.to_csv(output_het_file, index = False)

#Run the function
vcf_to_psmcfa(vcf_file, output_fq_file, output_het_file)
	
