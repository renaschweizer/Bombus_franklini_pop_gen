#!/bin/bash

#SBATCH -p atlas
#SBATCH --time=48:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=48   # 20 processor core(s) per node X 2 threads per core
#SBATCH --partition=atlas    # standard node(s)
#SBATCH --job-name="msprime"
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH -o msprime.%a.%A.out                    # File to which stdout will be written
#SBATCH -e msprime.%a.%A.err 
#SBATCH --account=beenome100

#This script simulates a demography WITH BOTTLENECK in msprime, then outputs a .vcf for the individuals
#specified in the "msprime_simulate_demography_bottleneck.py" script. A second python script,
#"make_consensus_diploid_sequence_from_vcf_one_chrom.py", is then called to generate consensus
#diploid sequence data for each individual, in the form of a fastq file. The invariant
#sequence data for the consensus diploid sequence is taken from chromosome one of the
#Bombus franklini reference genome.

#This shell script requires:
	#bcftools to be installed and in your PATH
	#msprime/tskit to be installed (conda environments are de/activated in the demography simulation section below)
	#Currently configured to work in the zsh shell

#make_consensus_diploid_sequence_from_vcf_one_chrom.py requires:
	#Path to (subset of) reference genome
	#pysam python module
	#Number of simulated chromosomes and chromosome length

#msprime_simulate_demography.py requires:
	#msprime python module
	#Mutation rate, recombination rate, chromosome length

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
# sbatch --array=n  ~/scripts/simulation_pipeline/msprime_simulate_demography_pipeline.slurm

module load  bcftools/1.19
module load miniconda3/24.3.0
source activate /project/beenome100/conda_envs/msprime

### Define some variables ###

wd=/90daydata/beenome100/rena_in_progress/bombus_franklini/msprime_outputs

#Number of individuals simulated
num_inds=25

#Name of demographic model being tested
model_name="bottleneck_model_26k_mut_rate_3.6e-9"

### Define some folder paths for results ###

#Define path to folder with python scripts
scripts_path="/home/rena.schweizer/scripts/simulation_pipeline/"

#Define file path for results and msprime-made vcf, which is then used for making individual vcf/fq files
results_folder=${wd}/results_${model_name}_strength_${strength}_${SLURM_ARRAY_TASK_ID}
mkdir -p ${results_folder}
msprime_vcf=${results_folder}/${model_name}.vcf

##############################################
##### Simulate the demography in msprime #####
##############################################

echo "Now running the demographic simulation in msprime"

#Run the msprime demographic simulation and generate a .vcf for multiple individuals
python ${scripts_path}/msprime_simulate_demography_one_chrom_bottleneck.py ${msprime_vcf} ${strength}

echo "Demographic simulation completed; the .vcf is at ${msprime_vcf}"

##########################################################
##### Convert msprime .vcf into individual .fq files #####
##########################################################

echo "Now converting the msprime .vcf to fastq inputs for f2psmcfa"

#Subset the vcf to get data for a single individual; can loop over multiple individuals here
#Create a list of the individuals that will have their own vcf and psmc files
#Make another variable for the for-loop, given that they are numbered with base-0
declare -i num_inds_loop 
num_inds_loop=${num_inds}-1

#A list to hold the names of the individuals to subset the vcf file
ind_list=()
for ((ind=0; ind<num_inds; ind++))
do
  ind_list+=(tsk_$ind)
done

#Make vcf and fastq PSMC-fa input files for each individual,
#then delete the .vcf files
het_table=${results_folder}/Heterozygosity_values_${model_name}_${strength}_${SLURM_ARRAY_TASK_ID}.csv
for indiv in "${ind_list[@]}"; do
	single_ind_vcf=${results_folder}/${indiv}.vcf
	single_ind_fq=${results_folder}/${indiv}.fq
	bcftools view -s ${indiv} ${msprime_vcf} > ${single_ind_vcf}
	python ${scripts_path}/make_consensus_diploid_sequence_from_vcf_one_chrom.py ${single_ind_vcf} ${single_ind_fq} ${het_table}
	gzip ${single_ind_fq}
	rm ${single_ind_vcf}
done

echo "Done converting .vcf files for f2psmcfa! They are in the ${results_folder} folder."
echo "The .vcf file made by msprime is at ${msprime_vcf}"
