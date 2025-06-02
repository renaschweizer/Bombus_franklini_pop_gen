#!/usr/bin/env python3
'''
This script simulates a specified demography in msprime and then writes out
a .vcf file for one chromosome. Specify the mutation rate, recombination rate,
chromosome length, and Ne of the simulation population at time 0 (current).'''

import math
import msprime
import os
import random
import sys

#Mutation rate for overlaying mutations on coalescent trees
mut_rate = 3.6e-9

#Recombination rate, per bp #8.7cM/Mb = 8.7e-8
#Having a high recombination rate slows down the simulations a lot
r_chrom = 8.7e-8

#How long is each chromosome
chr_length = int(1.3e7)

#Specify the path to the modified output vcf file
#E.g., python3 modify_vcf_to_ten_chroms.py one_ind_ten_chroms.vcf
output_file = sys.argv[1]

#Specify the path to the output file holding tree sequences
#output_tree_file = sys.argv[2]

#Random seed for tree sequence
ts_rand_seed = random.randint(1, 10000000)

#Give the effective pop. size of the current population (gen. 0)
current_haplodiploid_ne = 27500
#current_haplodiploid_ne = 100 # for testing

#Bottleneck strength, in generations of the coalescent
btl_str = int(sys.argv[2]) # 64000

#Correct for empirically estimated Ne, given haplodiploidy
current_diploid_ne = current_haplodiploid_ne * 1.33

#Set up the demography we want
#For the first 6 generations, the population size shrinks by 55% per generation (45% survival rate),
#so the coalescence probability gets higher as N gets smaller, where coal. prob. = 1 / 2N.
#Note that the coal. prob. changes, but the pop. size doesn't
demographic_model = msprime.Demography()
demographic_events=[
	demographic_model.add_population(initial_size = current_diploid_ne, growth_rate = 0.000012), #Size at t=0, then shrinks to ~11k individuals by gen. 40k
	demographic_model.add_instantaneous_bottleneck(time = 6, population = 0, strength = btl_str), #coalescence probability for pop. size 16459 ( = current_diploid_ne * 0.45)
	demographic_model.add_population_parameters_change(time=80000, growth_rate = -0.000038, population = 0), #Pop. goes from 11k - 110k by gen. 200k
	demographic_model.add_population_parameters_change(time=200000, growth_rate = 0, population = 0), #Stop growth at 200k gens ago
]


#Print the demographic events to the console to make sure the model is correct
dem_debug = demographic_model.debug()
print("The demography you simulated looks like:\n")
print(dem_debug)

#Set up the individuals to sample
#Sample 20 individuals 10 gens ago, and 5 individuals 60 gens ago
temporal_samples=[
	msprime.SampleSet(num_samples=6, population=0, time=0), #In the output vcf, these individuals should come first
	msprime.SampleSet(num_samples=1, population=0, time=3),
	msprime.SampleSet(num_samples=6, population=0, time=8),
	msprime.SampleSet(num_samples=1, population=0, time=18),
	msprime.SampleSet(num_samples=3, population=0, time=30),
	msprime.SampleSet(num_samples=5, population=0, time=31),
	msprime.SampleSet(num_samples=3, population=0, time=48),
]

#Run the actual demographic model and generate coalescent trees
ts = msprime.sim_ancestry(
	samples = temporal_samples,
	demography = demographic_model,
	recombination_rate = r_chrom,
	sequence_length = chr_length,
    random_seed = ts_rand_seed)

print("The coalescent trees you simulated look like:\n")
print(ts)

#Save the tree sequences to file
#ts.dump(output_tree_file)

#Mutate the tree sequence, and write to the vcf.
mts_rand_seed = random.randint(1, 10000000) #Random seed for adding mutations to the tree

#Add mutations to the chromosome
mts = msprime.sim_mutations(ts, rate = mut_rate, random_seed = mts_rand_seed)
    
#Calculate sample heterozygosity (pi, Nei and Li 1979) for the chromosome
print("Chromosome heterozygosity is", mts.diversity())

#Write out the variants to a .vcf file
with open(output_file, "w") as vcf_file:
	mts.write_vcf(vcf_file, allow_position_zero = True)
		

### OTHER USEFUL CODE BITS ###
'''#To see the variants and make sure they're there
for variant in mts.variants():
    print(
        "Variable site", variant.site.id,
        "at genome position", variant.site.position,
        ":", [variant.alleles[g] for g in variant.genotypes],
    )
    
#Code for plotting in R to see the pop. size change
#due to initial pop. size and growth rate
init_size=30000 #size at time 0, current population
a=0.000035 #growth rate
t=seq(0:30000) # num. gens. in the past for the population to shrink (grow forward in time)

pop_size = init_size*exp(-a*t)

plot(pop_size,
     ylab = "Pop. size",
     xlab = "Num. gens. ago")
'''
