#!/bin/bash

#SBATCH --time=48:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=48   
#SBATCH --job-name="related"
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH -o related.%a.%A.out                    

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
# usage: sbatch --export=species='bombus_franklini',vcf='vcf_name',prune='0.3' [or similar] ~/scripts/plink_relatedness.slurm

module load apptainer/1.1.9
module load plink2/2.00a4.3

wd=/90daydata/beenome100/rena_in_progress/${species}

ref=/project/beenome100/nuc_genomes_NCBI_reference/GCF_024516045.1_iyBomAffi1.2_genomic.fna

version='plink --version'

echo "calculate relatedness with plink ${version} on $(date)"

VCF_path=${wd}/filt_all/${vcf}.recode.vcf

# create list of pruned snps

plink2 --vcf ${VCF_path} --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 10 ${prune} --bad-ld --out ${wd}/relatedness/${vcf}
_${prune}

# calc relatedness using KING method, https://www.cog-genomics.org/plink/2.0/distance

plink2 --vcf ${VCF_path} --double-id --allow-extra-chr --set-missing-var-ids @:# --extract ${wd}/relatedness/${vcf}_${prune}.prune.in  --make-king-table
 --out ${wd}/relatedness/${vcf}_${prune}

mv related.* ${wd}/relatedness/std_outs/
