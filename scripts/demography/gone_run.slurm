#!/bin/bash

#SBATCH --time=48:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=48   # 20 processor core(s) per node X 2 threads per core
#SBATCH --job-name="gone_run"
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH -o gone_run.%a.%A.out                    # File to which stdout will be written
#SBATCH -e gone_run.%a.%A.err 

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
# usage: sbatch --export=species='bombus_franklini',vcf='vcf_name' ~/scripts/gone_run.slurm

wd=/90daydata/beenome100/rena_in_progress/${species}/gone_demography
cd ${wd}

echo "running gone on file ${vcf} on $(date)"
 
# run gone using gone file

bash script_GONE.sh ${vcf}

mkdir -p ${wd}/std_outs/
mv gone_run.* ${wd}/std_outs/
