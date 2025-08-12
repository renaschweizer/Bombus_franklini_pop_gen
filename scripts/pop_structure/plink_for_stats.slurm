#!/bin/bash

#SBATCH --time=48:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=48   # 20 processor core(s) per node X 2 threads per core
#SBATCH --job-name="stats"
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH -o stats.%a.%A.out                    # File to which stdout will be written
#SBATCH -e stats.%a.%A.err 

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
# usage: sbatch --export=species='bombus_franklini',vcf='vcf_name',type="gzvcf" ~/scripts/plink_for_stats.slurm
# this script will run plink on the designated file to generate a variety of statistical metrics


module load apptainer/1.1.9
module load plink2/2.00a4.3

wd=/90daydata/beenome100/rena_in_progress/${species}

ref=/project/beenome100/nuc_genomes_NCBI_reference/GCF_024516045.1_iyBomAffi1.2_genomic.fna

# Specify the path to the config file
config=${wd}/config_files/stats_to_run.txt # tab-delim file with array number and vcf stats, e.g. 1  idepth

# Extract the file name for the current $SLURM_ARRAY_TASK_ID
stat=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)

version='plink2 --version'

echo "calculate ${stat} with plink ${version}  on $(date)"

if [ $type == "gzvcf" ]
then
	vcf_path=${wd}/filt_all/$vcf".vcf.gz"
else
	vcf_path=${wd}/filt_all/$vcf".recode.vcf"
fi

# generate stat

plink2 --vcf ${vcf_path} --double-id --allow-extra-chr --set-missing-var-ids @:# --${stat} --out ${wd}/plink/${vcf}

mv stats.* ${wd}/plink/std_outs/
