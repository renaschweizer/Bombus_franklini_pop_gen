#!/bin/bash

#SBATCH --time=48:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=48   # 20 processor core(s) per node X 2 threads per core
#SBATCH --job-name="exractCOI"
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --output="std-%j-%N.out" # job standard output file (%j replaced by job id)
#SBATCH --error="std-%j-%N.err" # job standard error file (%j replaced by job id)

## THIS SCRIPT EXTRACTS COI SEQUENCE FROM AN ASSEMBLED GENOME IN FASTA FORMAT

## GET SET UP WITH CONDA
module load miniconda3/24.3.0

source activate /project/beenome100/conda_envs/phyluce-1.7.3

## RUN PHYLUCE TO SLICE OUT MATCHING SEQUENCE

phyluce_assembly_match_contigs_to_barcodes \
	--contigs ${my_fasta} \
	--barcodes ${my_barcodes} \
	--output ${outdir} \
	--no-bold

## SHUT DOWN
conda deactivate
