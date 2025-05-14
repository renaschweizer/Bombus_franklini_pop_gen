#!/bin/bash

#SBATCH -p atlas
#SBATCH --time=48:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=48   # 20 processor core(s) per node X 2 threads per core
#SBATCH --partition=bigmem    # standard node(s)
#SBATCH --job-name="qualimap"
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH -o qualimap.%a.%A.out                    # File to which stdout will be written
#SBATCH -e qualimap.%a.%A.err 
#SBATCH --account=beenome100

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

module load miniconda/4.12.0

source activate /project/beenome100/conda_envs/qualimap
export DISPLAY=:0 # turn off DISPLAY environment variable so that X11 system is not tested
export JAVA_OPTS="-Xmx1400g" 

wd=/90daydata/beenome100/rena_in_progress/${species}

# Specify the path to the config file
config=${wd}/config_files/trim_array.txt

# Extract the file name for the current $SLURM_ARRAY_TASK_ID
sample=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)

# Print to a file a message that includes the current $SLURM_ARRAY_TASK_ID, and the sample name
echo "This is array task ${SLURM_ARRAY_TASK_ID}, the file name is ${sample}."

version='qualimap --version'

echo "Assessing mapping quality via qualimap on $(date) with ${version}"

qualimap bamqc -nt 48 \
 -bam ${wd}/bam/${sample}_aln-pe.s.u.q.bam \
 -outdir ${wd}/mapping_quality_metrics/${sample}

mv qualimap.* ${wd}/mapping_quality_metrics/std_outs/

