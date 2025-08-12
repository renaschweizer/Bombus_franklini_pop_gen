#!/bin/bash

#SBATCH --mem-per-cpu=0
#SBATCH --time=168:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=48   # 20 processor core(s) per node X 2 threads per core
#SBATCH --job-name="spades"
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH -o spades.%a.%A.out                    # File to which stdout will be written
#SBATCH -e spades.%a.%A.err 

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

module load miniconda3/24.3.0

source activate /project/beenome100/conda_envs/spades

wd=/90daydata/beenome100/rena_in_progress/${species}

#### assembly with spades

echo "starting spades assembly for ${sample} on $(date)"

spades.py -1 ${wd}/spades_assembly/Bombus_franklini_${sample}_R1.trimmed.norm.fastq.gz \
	-2 ${wd}/spades_assembly/Bombus_franklini_${sample}_R2.trimmed.norm.fastq.gz \
	-o ${wd}/spades_assembly/${sample}_assembly \
	--threads 48

# if restarting
#spades.py -o ${wd}/spades_assembly \
#        --continue

echo "finishing spades assembly on $(date)"
