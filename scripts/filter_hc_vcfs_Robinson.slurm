#!/bin/bash

#SBATCH -p atlas
#SBATCH --time=48:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=48   # 20 processor core(s) per node X 2 threads per core
#SBATCH --partition=atlas    # standard node(s)
#SBATCH --job-name="filtVCFs"
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH -o filtJAR.%a.%A.out                    # File to which stdout will be written
#SBATCH -e filtJAR.%a.%A.err 
#SBATCH --account=beenome100

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

module load apptainer/1.1.6
module load singularity/3.8.3
module load singularity/ce-3.11.0
module load gatk4/4.1.8.0--py38h37ae868_0
module load vcftools/v0.1.16-1-deb_cv1

wd=/90daydata/beenome100/rena_in_progress/${species}

ref=/project/beenome100/nuc_genomes_NCBI_reference/GCF_024516045.1_iyBomAffi1.2_genomic.fna

# Specify the path to the config file
config=${wd}/config_files/genome_contigs_array.txt

# Extract the file name for the current $SLURM_ARRAY_TASK_ID
contig=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)

echo "Run filterVCFs script from Robinson et al. 2022 Science on  $(date)"

python3 ~/scripts/filterVCF_Robinson_et_al.py ${wd}/called_genotypes_HC/${species}_HC_geno_unfiltered_${contig}.vcf.gz | bgzip > ${wd}/called_genotypes_HC/${species}_HC_geno_annotatedJAR_${contig}.vcf.gz

echo "Finished annotating vcf file on $(date)"

echo "Filtering out sites that PASS or only have a WARNING on $(date)"

zcat ${wd}/called_genotypes_HC/${species}_HC_geno_annotatedJAR_${contig}.vcf.gz | awk '{ if ( $0~/^#/ || $7=="PASS" || $7=="WARN_missing"){print $0}}' \
| bgzip > ${wd}/filt_all/${species}_HC_geno_PASS_${contig}.vcf.gz

echo "Filtering out variant sites on $(date)"

zcat ${wd}/filt_all/${species}_HC_geno_PASS_${contig}.vcf.gz | awk '{ if ( $0~/^#/ || $8~"VariantType=SNP" ){print $0}}' \
| bgzip > ${wd}/filt_all/${species}_HC_geno_PASS_variants_${contig}.vcf.gz

mv filtJAR.* ${wd}/called_genotypes_HC/std_outs/

