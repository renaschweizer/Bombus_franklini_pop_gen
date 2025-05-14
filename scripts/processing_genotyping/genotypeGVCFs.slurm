#!/bin/bash

#SBATCH -p atlas
#SBATCH --time=48:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=48   # 20 processor core(s) per node X 2 threads per core
#SBATCH --partition=atlas    # standard node(s)
#SBATCH --job-name="genotypeGVCFs"
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH -o genotype.%a.%A.out                    # File to which stdout will be written
#SBATCH -e genotype.%a.%A.err 
#SBATCH --account=beenome100

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

module load apptainer/1.1.6
module load singularity/3.8.3
module load singularity/ce-3.11.0
module load gatk4/4.1.8.0--py38h37ae868_0

wd=/90daydata/beenome100/rena_in_progress/${species}

ref=/project/beenome100/nuc_genomes_NCBI_reference/GCF_024516045.1_iyBomAffi1.2_genomic.fna

# Specify the path to the config file
config=${wd}/config_files/genome_contigs_array.txt

# Extract the file name for the current $SLURM_ARRAY_TASK_ID
contig=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)

version='gatk --version'

echo "Run GATK GenotypeGVCFs ${version} on  $(date)"

gatk GenotypeGVCFs \
-R ${ref} \
-V gendb://${wd}/called_genotypes_HC/${species}_database_${contig} \
--intervals ${contig} \
--include-non-variant-sites \
-O ${wd}/called_genotypes_HC/${species}_HC_geno_unfiltered_${contig}.vcf.gz

mv genotype.* ${wd}/called_genotypes_HC/std_outs/

