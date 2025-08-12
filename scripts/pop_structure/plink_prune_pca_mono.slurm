#!/bin/bash

#SBATCH --time=48:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=48   # 20 processor core(s) per node X 2 threads per core
#SBATCH --job-name="pca"
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH -o pca.%a.%A.out                    # File to which stdout will be written
#SBATCH -e pca.%a.%A.err 

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
# usage: sbatch --export=species='bombus_franklini',vcf='vcf_name',prune='r2 value' ~/scripts/plink_prune_pca.slurm

module load apptainer/1.1.9
module load plink2/2.00a4.3

wd=/90daydata/beenome100/rena_in_progress/${species}

ref=/project/beenome100/nuc_genomes_NCBI_reference/GCF_024516045.1_iyBomAffi1.2_genomic.fna

version='plink2 --version'

echo "generate pca with plink ${version} and prune parameter ${prune} on $(date)"

VCF_path=${wd}/filt_all/${vcf}.recode.vcf
 
# create list of pruned snps

plink2 --vcf ${VCF_path} --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 10 ${prune} --bad-ld --out ${wd}/plink/${vcf}_${prune}

# generate frequency data

plink2 --vcf ${VCF_path} --double-id --allow-extra-chr --set-missing-var-ids @:# --freq --out ${wd}/plink/${vcf}_${prune}

# run pca

plink2 --vcf ${VCF_path} --double-id --allow-extra-chr --set-missing-var-ids @:# --extract ${wd}/plink/${vcf}_${prune}.prune.in --make-bed --pca --out ${wd}/plink/${vcf}_${prune} --read-freq ${wd}/plin
k/${vcf}_${prune}.afreq


mv pca.* ${wd}/plink/std_outs/
