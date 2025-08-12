#!/bin/bash

#SBATCH --time=48:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=48   # 20 processor core(s) per node X 2 threads per core
#SBATCH --job-name="split"
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH -o split.%a.%A.out                    # File to which stdout will be written
#SBATCH -e split.%a.%A.err 

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

# use psmc tools to split input file 

version='${software}/psmc/psmc --version'

echo "Split psmc input file for ${sample} on $(date) with ${version}"

i=${wd}/psmc/${dir}/${sample}.diploid.psmcfa

${software}/psmc/utils/splitfa ${i} > ${wd}/psmc/${dir}/${sample%.psmcfa}_split.psmcfa

echo "done splitting psmc file on $(date)"

mv split.* ${wd}/psmc/${dir}/std_outs/
