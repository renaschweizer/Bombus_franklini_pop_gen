#!/bin/bash

# usage: sbatch --array=19-66 --export=species='bombus_affinis' ~/scripts/bwa_array.slurm
# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

module load apptainer/1.1.6
module load singularity/3.8.3
module load singularity/ce-3.11.0
module load bwa/0.7.17
module load samtools/v1.9-4-deb_cv1

wd=/90daydata/beenome100/rena_in_progress/${species}

ref=/project/beenome100/nuc_genomes_NCBI_reference/GCF_024516045.1_iyBomAffi1.2_genomic.fna

# Specify the path to the config file
config=${wd}/config_files/trim_array.txt

# Extract the file name for the current $SLURM_ARRAY_TASK_ID
sample=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)

# Print to a file a message that includes the current $SLURM_ARRAY_TASK_ID, and the sample name
echo "This is array task ${SLURM_ARRAY_TASK_ID}, the file name is ${sample}."

 version='bwa --version'

echo "Generating mapped reads on $(date) with ${version}"

bwa mem -M -t 48 \
${ref} \
${wd}/trimmed_fastq/${sample}_R1.trimmed.fastq \
${wd}/trimmed_fastq/${sample}_R2.trimmed.fastq \
> ${wd}/bam/${sample}_aln-pe.sam

# convert SAM to BAM

 version='samtools --version'

echo "Generating BAM file from SAM on $(date) with ${version}"

samtools view -b -h -S ${wd}/bam/${sample}_aln-pe.sam \
	> ${wd}/bam/${sample}_aln-pe.bam

echo "Fix mates and fill in insert sizes, then sort on $(date)" 

samtools fixmate -m \
${wd}/bam/${sample}_aln-pe.bam \
${wd}/bam/${sample}_aln-pe.fix.bam

samtools sort  -@ 48 \
-o ${wd}/bam/${sample}_aln-pe.s.bam \
${wd}/bam/${sample}_aln-pe.fix.bam

echo "Remove PCR duplicates on $(date)" 

samtools markdup -r ${wd}/bam/${sample}_aln-pe.s.bam \
${wd}/bam/${sample}_aln-pe.s.u.bam

echo "Remove low quality reads on $(date)"

samtools index -@ 48  ${wd}/bam/${sample}_aln-pe.s.u.bam 

samtools view -h -b -q 20 ${wd}/bam/${sample}_aln-pe.s.u.bam \
> ${wd}/bam/${sample}_aln-pe.s.u.q.bam


mv bwa.* ${wd}/bam/std_outs/

