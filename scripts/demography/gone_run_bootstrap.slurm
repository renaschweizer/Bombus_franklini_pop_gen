#!/bin/bash

#SBATCH --time=48:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=48   # 20 processor core(s) per node X 2 threads per core
#SBATCH --job-name="gone_run"
#SBATCH --mail-type=FAIL
#SBATCH -o gone_run.%a.%A.out                    # File to which stdout will be written
#SBATCH -e gone_run.%a.%A.err 

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
# usage: sbatch --export=species='bombus_franklini',vcf='vcf_name',boot=1 ~/scripts/gone_run.slurm

wd=/90daydata/beenome100/rena_in_progress/${species}/gone_demography/bootstrap_${boot}

echo "running gone on file ${vcf} for bootstrap ${boot} on $(date)"
 
# run gone using gone file

cd ${wd}

cp -r /project/beenome100/software/GONE/Linux/PROGRAMMES/ . 
cp -r /project/beenome100/software/GONE/Linux/script_GONE.sh . 
cp /90daydata/beenome100/rena_in_progress/${species}/gone_demography/INPUT_PARAMETERS_FILE . 

chmod a+rwx */*

bash script_GONE.sh ${vcf}_${boot}

mkdir -p ../std_outs/
mv ../gone_run.* ../std_outs/
