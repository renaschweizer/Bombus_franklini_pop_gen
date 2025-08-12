#!/bin/bash

#SBATCH --time=48:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=48   # 20 processor core(s) per node X 2 threads per core
#SBATCH --job-name="counts"
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH -o count.%a.%A.out                    # File to which stdout will be written
#SBATCH -e count.%a.%A.err 

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
# usage: sbatch --export=species='bombus_franklini',vcf='vcf_name',type="gzvcf" ~/scripts/plink_counts.slurm

module load apptainer/1.1.9
module load plink2/2.00a4.3

wd=/90daydata/beenome100/rena_in_progress/${species}

ref=/project/beenome100/nuc_genomes_NCBI_reference/GCF_024516045.1_iyBomAffi1.2_genomic.fna

version='plink2 --version'

echo "calculate counts of homozygous and heterozygous sites with plink ${version}  on $(date)"

if [ $type == "gzvcf" ]
then
	vcf_path=${wd}/filt_all/$vcf".vcf.gz"
else
	vcf_path=${wd}/filt_all/$vcf".recode.vcf"
fi

# generate count data

plink2 --vcf ${vcf_path} --double-id --allow-extra-chr --set-missing-var-ids @:# --sample-counts cols=hom,het --out ${wd}/plink/${vcf}_het

mv count.* ${wd}/plink/std_outs/
