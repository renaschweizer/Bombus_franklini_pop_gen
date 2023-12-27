#!/bin/bash

#SBATCH -p atlas
#SBATCH --time=48:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=48   # 20 processor core(s) per node X 2 threads per core
#SBATCH --partition=atlas    # standard node(s)
#SBATCH --job-name="map.stats"
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH -o stats.%a.%A.out                    # File to which stdout will be written
#SBATCH -e stats.%a.%A.err 
#SBATCH --account=beenome100

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

module load apptainer/1.1.6
module load singularity/3.8.3
module load singularity/ce-3.11.0
module load samtools/v1.9-4-deb_cv1

wd=/90daydata/beenome100/rena_in_progress/${species}

# Specify the path to the config file
config=${wd}/config_files/trim_array.txt

# Extract the file name for the current $SLURM_ARRAY_TASK_ID
sample=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)

# Print to a file a message that includes the current $SLURM_ARRAY_TASK_ID, and the sample name
echo "This is array task ${SLURM_ARRAY_TASK_ID}, the file name is ${sample}."

echo "Counting num raw and filtered reads on $(date)" 
echo "sampleID	rawR1	rawR2	trimR1	trimR2	mapPR	fixmatePR	uniqPR	MQ20PR" > ${wd}/mapping_quality_metrics/${sample}_mapping_stats.txt
echo ${sample} >> ${wd}/mapping_quality_metrics/${sample}_mapping_stats.txt
echo $(zcat ${wd}/raw_data/${sample}_R1.fastq.gz|wc -l)/4|bc >> ${wd}/mapping_quality_metrics/${sample}_mapping_stats.txt
echo $(zcat ${wd}/raw_data/${sample}_R2.fastq.gz|wc -l)/4|bc >> ${wd}/mapping_quality_metrics/${sample}_mapping_stats.txt

echo $(cat ${wd}/trimmed_fastq/${sample}_R1.trimmed.fastq| wc -l)/4|bc >> ${wd}/mapping_quality_metrics/${sample}_mapping_stats.txt
echo $(cat ${wd}/trimmed_fastq/${sample}_R2.trimmed.fastq| wc -l)/4|bc >> ${wd}/mapping_quality_metrics/${sample}_mapping_stats.txt

 version='samtools --version'

echo "Counting number of mapped, fix mate pair, and non-dup reads via samtools on $(date) with ${version}"

samtools view -c -f 0x2 ${wd}/bam/${sample}_aln-pe.bam >> ${wd}/mapping_quality_metrics/${sample}_mapping_stats.txt
samtools view -c -f 0x2 ${wd}/bam/${sample}_aln-pe.fix.bam >>  ${wd}/mapping_quality_metrics/${sample}_mapping_stats.txt
samtools view -c -f 0x2 ${wd}/bam/${sample}_aln-pe.s.u.bam >>  ${wd}/mapping_quality_metrics/${sample}_mapping_stats.txt
samtools view -c -f 0x2 ${wd}/bam/${sample}_aln-pe.s.u.q.bam >> ${wd}/mapping_quality_metrics/${sample}_mapping_stats.txt

tr '\n' '\t' <${wd}/mapping_quality_metrics/${sample}_mapping_stats.txt | sed 's/Baff/\nBaff/g' | sed 's/Bombus/\nBombus/g' > ${wd}/mapping_quality_metrics/${sample}_tmp.txt
mv ${wd}/mapping_quality_metrics/${sample}_tmp.txt ${wd}/mapping_quality_metrics/${sample}_mapping_stats.txt

mv stats.* ${wd}/mapping_quality_metrics/std_outs/

