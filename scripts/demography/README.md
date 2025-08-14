Scripts to run demography analyses

Following recommendations of PSMC [github](https://github.com/lh3/psmc) : 

+ `psmc.slurm` to generate input files and run PSMC for an individual

+ `psmc_splitfa.slurm` and `psmc_bootstrap.slurm` to split input files and run multiple bootsraps of PSMC

Following recommendations of GONE [github](https://github.com/esrud/GONE)

+ `gone_prep_input.slurm` and `gone_run.slurm` to run GONE for all individuals
+ `gone_prep_input_bootstrap.slurm` and `gone_run_bootstrap.slurm` to randomly subsample to 15 individuals and run GONE for iterative bootstraps
