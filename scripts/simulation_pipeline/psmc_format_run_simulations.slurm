#!/bin/bash

#SBATCH -p atlas
#SBATCH --time=48:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=48   # 20 processor core(s) per node X 2 threads per core
#SBATCH --partition=atlas    # standard node(s)
#SBATCH --job-name="psmc"
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH -o psmc.%a.%A.out                    # File to which stdout will be written
#SBATCH -e psmc.%a.%A.err 
#SBATCH --account=beenome100

# usage: sbatch --array=19-66 --export=species='bombus_affinis' ~/scripts/psmc_format_run.slurm
# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

module load miniconda3/24.3.0

#wd=/90daydata/beenome100/rena_in_progress/${species}
conda_envs=/project/beenome100/conda_envs/
software=/project/beenome100/software/

ref=/project/beenome100/nuc_genomes_NCBI_reference/GCF_024516045.1_iyBomAffi1.2_genomic.fna

# activate psmc environment
source activate ${conda_envs}/psmc

# Specify the path to the config file
#config=${wd}/config_files/depth_filter_calcs.txt

# Extract the file name for the current $SLURM_ARRAY_TASK_ID
#sample=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)

# Print to a file a message that includes the current $SLURM_ARRAY_TASK_ID, and the sample name
#echo "This is array task ${SLURM_ARRAY_TASK_ID}, the file name is ${sample}."

# use psmc tools to generate input file from consensus sequence

version='${software}/psmc/psmc --version'

#echo "Generating psmc input file for ${sample} on $(date) with ${version}"

#${software}/psmc/utils/fq2psmcfa -q20  ${wd}/psmc/${sample}.diploid.fq.gz > ${wd}/psmc/${sample}.diploid.psmcfa

#echo "done reformatting input file on $(date)" 

#echo "running psmc on ${sample} on $(date)"

${software}/psmc/psmc -N30 -t15 -r5 -p "4+20*2+6" -o ten_chroms.psmc_output.psmc ten_chroms_psmcfa_input.txt
echo "done running psmc on $(date)"

#mv psmc.* ${wd}/psmc/std_outs/

