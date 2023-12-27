These are the scripts used to process the raw sequence data, map to the reference genome, call genotypes, and filter genotypes. 

#### Make config files

e.g., raw_fastq_array.txt

```{unix}
ArrayTaskID   FileName
1 raw_data/Baff_256/Baff_256_R1.fastq.gz
2 raw_data/Baff_256/Baff_256_R2.fastq.gz
3 raw_data/Baff_257/Baff_257_R1.fastq.gz
4 raw_data/Baff_257/Baff_257_R2.fastq.gz
```

e.g., trim_array.txt
```{unix}
ArrayTaskID   FileName
1 Baff_256
2 Baff_257
3 Baff_258
4 Baff_ESI3Aug211
```

e.g., fastq_trimmed_array.txt
```{unix}
ArrayTaskID   FileName
1 trimmed_fastq/Baff_256_R1.trimmed.fastq
2 trimmed_fastq/Baff_256_R2.trimmed.fastq
3 trimmed_fastq/Baff_257_R1.trimmed.fastq
4 trimmed_fastq/Baff_257_R2.trimmed.fastq
```

e.g., addRG_info_array.txt
```{unix}
ArrayTaskID	FileName	LaneNum	RunID
1	Baff_256	1	707
2	Baff_257	1	707
3	Baff_258	1	707
4	Baff_ESI3Aug211	1	707
```

e.g., genome_contigs_array.txt
```{unix}
ArrayTaskID   FileName
1 NC_066344.1
2 NC_066345.1
3 NC_066346.1
4 NC_066347.1
```


#### Run fastqc on raw data. 

```{unix}
sbatch --array=37-132 --export=species='bombus_affinis' ~/scripts/fastqc_array.slurm
```

#### Trim adapters. 
```{unix}
sbatch --array=19-66 --export=species='bombus_affinis' ~/scripts/trimmomatic_array.slurm
```

#### Map trimmed reads to reference genome.
```{unix}
sbatch --array=19-66 --export=species='bombus_affinis' ~/scripts/bwa_array.slurm 
```

#### Run mapping stats.
```{unix}
 sbatch --array=19-66 --export=species='bombus_affinis'  ~/scripts/mapping_stats.slurm
 sbatch --array=19-66 --export=species='bombus_affinis' ~/scripts/qualimap.slurm 
```
 
#### Double check quality of trimmed fastq reads
```{unix}
sbatch --array=37-132 --export=species='bombus_affinis' ~/scripts/fastqc_array.slurm
```

#### AddRG information. 
```{unix}
sbatch --array=19-66 --export=species='bombus_affinis' ~/scripts/addRG_array.slurm
```

#### Run BQSR from individual bams. I've split up the runs by individual contigs (n=858 total).
```{unix}
sbatch --array=1-858 --export=species='bombus_affinis' ~/scripts/bqsr_from_ibams_array.slurm
```

#### Gather/Apply BQSR results.
```{unix}
sbatch --export=round=1,species='bombus_affinis' ~/scripts/gather_applybqsr.slurm
```

#### Round 2 of BQSR.
```{unix}
sbatch --array=1-858 --export=species='bombus_affinis' ~/scripts/bqsr_round2_array.slurm
```

#### Gather/apply from round 2. 
```{unix}
sbatch --export=round=2,species='bombus_affinis' ~/scripts/gather_applybqsr.slurm
```

#### Call HaplotypeCaller on each individual.
```{unix}
sbatch --array=19-66 --export=species='bombus_affinis' ~/scripts/haplotype_caller.slurm
```

#### Run GenomicsDBImport in array across contigs.
```{unix}
sbatch --array=1-858 --export=species='bombus_affinis' ~/scripts/genomicsDBimport_array.slurm
```

#### Genotype GVCFs
```{unix}
sbatch --array=1-858 --export=species='bombus_affinis' ~/scripts/genotypeGVCFs.slurm
```

#### Merge VCF files.
```{unix}
ls filt_SNPs/*.recode.vcf | awk '{print "/90daydata/beenome100/rena_in_progress/bombus_affinis/"$1}' > filt_SNPs/filt_snp_files.txt
sbatch --export=species='bombus_affinis' ~/scripts/merge_vcfs.slurm

ls filt_all/*.recode.vcf | awk '{print "/90daydata/beenome100/rena_in_progress/bombus_affinis/"$1}' > filt_all/filt_all_files.txt
sbatch --export=species='bombus_affinis' ~/scripts/merge_vcfs_mono.slurm
```


