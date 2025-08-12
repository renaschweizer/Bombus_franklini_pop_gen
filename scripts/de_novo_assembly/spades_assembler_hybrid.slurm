#!/bin/bash

#SBATCH --time=168:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=48   
#SBATCH --mem-per-cpu=0
#SBATCH --job-name="spades"
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH -o spades.%a.%A.out                    # File to which stdout will be written
#SBATCH -e spades.%a.%A.err 

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

module load miniconda3/24.3.0

source activate /project/beenome100/conda_envs/spades

wd=/90daydata/beenome100/rena_in_progress/bombus_franklini/spades_assembly

#### assembly with spades

echo "starting spades assembly on $(date)"

spades.py -1 ${wd}/Bombus_franklini_GNS104_R1.trimmed.norm.fastq \
	-2 ${wd}/Bombus_franklini_GNS104_R2.trimmed.norm.fastq \
	-o ${wd}/hybrid \
	--threads 48 \
	--nanopore /90daydata/beenome100/rena_in_progress/bombus_franklini/raw_data/nanopore/Bombus_franklini_combined_simplex_concat.fastq.gz 

# if restarting
#spades.py -o ${wd}/spades_assembly \
#        --continue

echo "finishing spades assembly on $(date)"
