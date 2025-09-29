# RoyerTUAssignments
import allel
import pandas as pd
import numpy as np

tmpvcf = 'C:/Users/laure/RoyerTUAssignments/HumanEvolutionPt2/test_500.vcf.gz'
longvcf = 'C:/Users/laure/Documents/BIOL5131/Week4/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz' ##update to just 1000lines to run from GitHub repo

### read in pop data
panel = pd.read_csv('integrated_call_samples_v3.20130502.ALL.panel', sep='\t')
### chop it up
pop_filter = None
super_pop_filter = 'EUR'

# select the sample names from metadata
if pop_filter is not None:
    samples_to_keep = panel.loc[panel['pop'] == pop_filter, 'sample'].tolist()
elif super_pop_filter is not None:
    samples_to_keep = panel.loc[panel['super_pop'] == super_pop_filter, 'sample'].tolist()
else:
    samples_to_keep = panel['sample'].tolist()  # use all samples


###adding in pop filters
callset = allel.read_vcf(tmpvcf, fields=['calldata/GT','variants/REF', 'variants/ALT', 'variants/POS'], samples=samples_to_keep)

genotypes = allel.GenotypeArray(callset['calldata/GT'])
alt_alleles = callset['variants/ALT'] ##pull from last week down
ref_alleles = callset['variants/REF']
allele_pos = callset['variants/POS']

###allele counts/SNP, replacing .is_snp and _is_biallel
allele_counts = genotypes.count_alleles() ##this does not skip >2 per, need to add manual filter as df
allele_counts = allel.AlleleCountsArray(allele_counts)

is_snp = np.array([len(r) == 1 and len(a[0]) == 1 for r, a in zip(ref_alleles, alt_alleles)])
is_biallelic = allele_counts.max_allele() ==1

fullfilter = is_snp & is_biallelic

genotypes = genotypes[fullfilter]
alt_alleles = alt_alleles[fullfilter]
ref_alleles = ref_alleles[fullfilter]
allele_pos = allele_pos[fullfilter]
allele_counts= allele_counts[fullfilter]

###bring here for filtered, delete up top if no diff
het_counts = genotypes.count_het(axis=1) ##count_het built in from allel, can rewrite to extract df and count per SNP
n_samples = genotypes.n_samples ##use len here when swapping
obs_het = het_counts / n_samples ##this part not automated, can plug here

allele_freqs = allele_counts.to_frequencies() ###convert to frequencies also manual after

## HW dom rec combo for mask
p = allele_freqs[:, 0] ##dom
q = allele_freqs[:, 1] ##rec, sum 1
exp_het = 2 * p * q

diff = obs_het - exp_het

results = pd.DataFrame({
    'observed_het': obs_het,
    'expected_het': exp_het,
    'difference': diff,
    'position': allele_pos,
    'ref': ref_alleles,
    'alt': [a[0] for a in alt_alleles]
})

print(results.head()) ##confirm <1
reportname = "pop-" + str(pop_filter) + "_super-" + str(super_pop_filter) + "_report.csv"
results.to_csv(reportname)