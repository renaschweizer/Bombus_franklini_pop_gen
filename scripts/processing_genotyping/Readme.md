These are the scripts used to process the raw sequence data, map to the reference genome, call genotypes, and filter genotypes on the SciNet ATLAS cluster. For more detailed information on the GATK workflow, please see [my other github repo on GATK](https://github.com/renaschweizer/congen-gatk).

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
R1 and R2 run separately so array=1-2n where n is # individuals sequenced, assuming one pair of reads per individual.

```{unix}
sbatch --array=37-132 --export=species='bombus_affinis' ~/scripts/fastqc_array.slurm
```

#### Trim adapters. 
Array is per individual. 

```{unix}
sbatch --array=19-66 --export=species='bombus_affinis' ~/scripts/trimmomatic_array.slurm
```

#### Map trimmed reads to reference genome.
Array is per individual. 

```{unix}
sbatch --array=19-66 --export=species='bombus_affinis' ~/scripts/bwa_array.slurm 
```

#### Run mapping stats.
Array is per individual. 

```{unix}
 sbatch --array=19-66 --export=species='bombus_affinis'  ~/scripts/mapping_stats.slurm
 sbatch --array=19-66 --export=species='bombus_affinis' ~/scripts/qualimap.slurm 
```
 
#### Double check quality of trimmed fastq reads
R1 and R2 run separately so array=1-2n where n is # individuals sequenced, assuming one pair of reads per individual.

```{unix}
sbatch --array=37-132 --export=species='bombus_affinis' ~/scripts/fastqc_array.slurm
```

#### AddRG information. 
Array is per individual.

```{unix}
sbatch --array=19-66 --export=species='bombus_affinis' ~/scripts/addRG_array.slurm
```

#### Run BQSR from individual bams. I've split up the runs by individual contigs (n=858 total).
Array is per genome contig.

```{unix}
sbatch --array=1-858 --export=species='bombus_affinis' ~/scripts/bqsr_from_ibams_array.slurm
```

#### Gather/Apply BQSR results.
No array.

```{unix}
sbatch --export=round=1,species='bombus_affinis' ~/scripts/gather_applybqsr.slurm
```

#### Round 2 of BQSR.
Array is per genome contig.

```{unix}
sbatch --array=1-858 --export=species='bombus_affinis' ~/scripts/bqsr_round2_array.slurm
```

#### Gather/apply from round 2. 
No array.

```{unix}
sbatch --export=round=2,species='bombus_affinis' ~/scripts/gather_applybqsr.slurm
```

#### Call HaplotypeCaller on each individual.
Array is per individual.

```{unix}
sbatch --array=19-66 --export=species='bombus_affinis' ~/scripts/haplotype_caller.slurm
```

#### Run GenomicsDBImport in array across contigs.
Array is per genome contig. 

```{unix}
sbatch --array=1-858 --export=species='bombus_affinis' ~/scripts/genomicsDBimport_array.slurm
```

#### Genotype GVCFs
Array is per genome contig. 

```{unix}
sbatch --array=1-858 --export=species='bombus_affinis' ~/scripts/genotypeGVCFs.slurm
```

#### Filter VCF files. 
Array is per genome contig. Follows GATK Best Practices for hard filtering, plus a min (1/3 mean) and max (2x mean) DP, excess heterozygosity, following Robinson, et al. 2022, Science.

```{unix}
sbatch --array=1-858 --export=species='bombus_affinis' ~/scripts/filter_hc_vcfs_Robinson.slurm
```

#### Merge VCF files.

```{unix}
ls bombus_affinis_HC_geno_PASS_N* > pass_files.list
 ls bombus_affinis_HC_geno_PASS_variants_N* > pass_variant_files.list

picard MergeVcfs I=pass_files.list O=bombus_affinis_combined_PASS.vcf.gz
picard MergeVcfs I=pass_variant_files.list O=bombus_affinis_combined_PASS_variants.vcf.gz
```



