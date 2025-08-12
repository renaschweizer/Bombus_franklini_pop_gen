#!/bin/bash

#SBATCH --time=48:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=48   # 20 processor core(s) per node X 2 threads per core
#SBATCH --job-name="boot"
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH -o boot.%a.%A.out                    # File to which stdout will be written
#SBATCH -e boot.%a.%A.err 
# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

module load miniconda3/24.3.0

wd=/90daydata/beenome100/rena_in_progress/${species}
conda_envs=/project/beenome100/conda_envs/
software=/project/beenome100/software/

# activate psmc environment
source activate ${conda_envs}/psmc

# Specify the path to the config file
config=${wd}/config_files/depth_filter_calcs.txt

# Extract the file name for the current $SLURM_ARRAY_TASK_ID
sample=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)

# Print to a file a message that includes the current $SLURM_ARRAY_TASK_ID, and the sample name
echo "This is array task ${SLURM_ARRAY_TASK_ID}, the file name is ${sample}."

# run psmc with bootsrapping 
# use -b flag 

version='${software}/psmc/psmc --version'

echo "Run psmc bootstrap ${iter} for ${sample} on $(date) with ${version}"

i=${wd}/psmc/${dir}/${sample}split.psmcfa

${software}/psmc/psmc -b -N30 -t15 -r5 -p "4+20*2+6" -o ${wd}/psmc/${dir}/${sample%.psmcfa}_split_boot_${iter}.out.psmc  ${wd}/psmc/${dir}/${sample%.psmcfa}_split.psmcfa

echo "done running psmc bootstrap file on $(date)"

mv boot.* ${wd}/psmc/${dir}/std_outs/
