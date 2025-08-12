#!/bin/bash

#SBATCH --time=48:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=48   # 20 processor core(s) per node X 2 threads per core
#SBATCH --job-name="admixture"
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH -o admix.%a.%A.out                    # File to which stdout will be written
#SBATCH -e admix.%a.%A.err 

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
# usage: sbatch --export=species='bombus_franklini',vcf='vcf_name',K='[n] ~/scripts/admixture.slurm
# K is the number of pops in admixture

module load apptainer/1.1.9
module load miniconda3/24.3.0

wd=/90daydata/beenome100/rena_in_progress/${species}

version='adxmiture --version'

echo " run admixture ${version} on $(date)"

source activate /project/beenome100/conda_envs/admixture

admixture ${wd}/plink/${vcf}_fixChr.bed ${K}

conda deactivate

mv admix.* ${wd}/admixture/std_outs/
