Scripts and files used to run simulation pipeline with and without recent pathogen bottleneck in msprime. 

Setup to run in an array format for simulation repetitions, e.g. 

```{unix}
sbatch --array=n  ~/scripts/simulation_pipeline/msprime_simulate_demography_pipeline.slurm
```

or 

```{unix}
sbatch --array=n  ~/scripts/simulation_pipeline/msprime_simulate_demography_pipeline_bottleneck.slurm
```

Both scripts call a separate python script that runs msprime, as well as a python script that creates psmc input files. 

