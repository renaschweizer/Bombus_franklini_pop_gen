#!/bin/bash

#SBATCH --time=48:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=48   
#SBATCH --job-name="mergeVCFs"
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH -o merge.%a.%A.out                    
#SBATCH -e merge.%a.%A.err 

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
# Combine vcf files
module load apptainer/1.1.6
module load singularity/3.8.3
module load singularity/ce-3.11.0
module load picard/2.23.2--0

wd=/90daydata/beenome100/rena_in_progress/${species}

ref=/project/beenome100/nuc_genomes_NCBI_reference/GCF_024516045.1_iyBomAffi1.2_genomic.fna

version='picard --version'

echo "Run picard MergeVcfs ${version} on  $(date)"

picard MergeVcfs \
I= ${wd}/filt_SNPs/filt_snp_files.txt \
O=${wd}/filt_SNPs/${species}_combined_filt.vcf.gz


mv merge.* ${wd}/filt_SNPs/std_outs/

