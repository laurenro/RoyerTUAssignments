# RoyerTUAssignments
import allel
import pandas as pd
# import collections
# import matplotlib.pyplot as plt
import numpy as np

tmpvcf = 'C:/Users/laure/RoyerTUAssignments/HumanEvolutionPt2/test_500.vcf.gz'
longvcf = 'C:/Users/laure/Documents/BIOL5131/Week4/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz' ##update to just 1000lines to run from GitHub repo

# # readvcf = allel.read_vcf(longvcf)
# # keys = readvcf.keys()
# # print(keys) 

###just checking to make sure VCF genotypes check out
callset = allel.read_vcf(tmpvcf, fields=["calldata/GT"])
genotypes = allel.GenotypeArray(callset["calldata/GT"])
# print(genotypes.shape) ##check for size of chopped vcf, can remove for submission

###use allel from scikit allal to get heterozygote counts per variant
het_counts = genotypes.count_het(axis=1) ##count_het built in from allel, can rewrite to extract df and count per SNP
n_samples = genotypes.n_samples ##use len here when swapping
obs_het = het_counts / n_samples ##this part not automated, can plug here

###allele counts/SNP
allele_counts = genotypes.count_alleles() ##this does not skip >2 per, need to add manual filter as df
allele_freqs = allele_counts.to_frequencies() ###convert to frequencies also manual after

## HW dom rec combo for mask
p = allele_freqs[:, 0] ##dom
q = allele_freqs[:, 1] ##rec, sum 1
exp_het = 2 * p * q

diff = obs_het - exp_het

results = pd.DataFrame({
    "observed_het": obs_het,
    "expected_het": exp_het,
    "difference": diff
})

print(results.head()) ##confirm <1
