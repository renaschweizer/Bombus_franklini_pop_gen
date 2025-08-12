#!/bin/bash

#SBATCH --time=48:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=48   # 20 processor core(s) per node X 2 threads per core
#SBATCH --job-name="gone_prep"
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH -o gone_prep.%a.%A.out                    # File to which stdout will be written
#SBATCH -e gone_prep.%a.%A.err 

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
# usage: sbatch --export=species='bombus_franklini',vcf='vcf_name' ~/scripts/gone_prep_input.slurm

module load plink2/2.00a4.3

wd=/90daydata/beenome100/rena_in_progress/${species}

version='plink2 --version'

echo "generate map input files with plink ${version} on $(date)"

VCF_path=${wd}/filt_all/${vcf}.recode.vcf
 
# create .ped and .map output with 18 chromosomes

plink2 --vcf ${VCF_path} --double-id --allow-extra-chr --set-missing-var-ids @:# --export ped --chr NC_066344.1,NC_066345.1,NC_066346.1,NC_066347.1,NC_066348.1,NC_066349.1,NC_066350.1,NC_066351.1,NC_06
6352.1,NC_066353.1,NC_066354.1,NC_066355.1,NC_066356.1,NC_066357.1,NC_066358.1,NC_066359.1,NC_066360.1,NC_066361.1  --out ${wd}/gone_demography/${vcf}

cat ${wd}/gone_demography/${vcf}.map | sed 's/NC_066344.1/1/g' | sed 's/NC_066345.1/2/g' | sed 's/NC_066346.1/3/g' | sed 's/NC_066347.1/4/g' | sed 's/NC_066348.1/5/g' | sed 's/NC_066349.1/6/g' | sed 's
/NC_066350.1/7/g' | sed 's/NC_066351.1/8/g' | sed 's/NC_066352.1/9/g' | sed 's/NC_066353.1/10/g' | sed 's/NC_066354.1/11/g' | sed 's/NC_066355.1/12/g' | sed 's/NC_066356.1/13/g' | sed 's/NC_066357.1/14
/g' | sed 's/NC_066358.1/15/g' | sed 's/NC_066359.1/16/g' | sed 's/NC_066360.1/17/g' | sed 's/NC_066361.1/18/g' > ${wd}/gone_demography/tmp.map

mv ${wd}/gone_demography/tmp.map ${wd}/gone_demography/${vcf}.map 


mkdir -p ${wd}/gone_demography/std_outs/
mv gone_prep.* ${wd}/gone_demography/std_outs/
