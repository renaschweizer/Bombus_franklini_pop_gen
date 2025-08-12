#!/bin/bash

# usage: sbatch --array=19-66 --export=species='bombus_affinis' ~/scripts/psmc_format_run.slurm
# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

module load miniconda3/24.3.0

conda_envs=/project/beenome100/conda_envs/
software=/project/beenome100/software/

ref=/project/beenome100/nuc_genomes_NCBI_reference/GCF_024516045.1_iyBomAffi1.2_genomic.fna

# activate psmc environment
source activate ${conda_envs}/psmc

${software}/psmc/psmc -N30 -t15 -r5 -p "4+20*2+6" -o ten_chroms.psmc_output.psmc ten_chroms_psmcfa_input.txt
echo "done running psmc on $(date)"
