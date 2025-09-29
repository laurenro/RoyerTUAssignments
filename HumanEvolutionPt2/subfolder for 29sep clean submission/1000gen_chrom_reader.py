# RoyerTUAssignments
### note for future this is build for single chrom only
import allel ##ref here: https://scikit-allel.readthedocs.io/en/stable/
import pandas as pd
import numpy as np

tmpvcf = 'test_500.vcf.gz'
# longvcf = 'ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz' ##cannot upload full file in github repo so defaulted to tmp

### read in pop data
panel = pd.read_csv('integrated_call_samples_v3.20130502.ALL.panel', sep='\t')
### chop it up
pop_filter = ['YRI','CEU','PJL']
# summary_list=[]
df_list=[]
# super_pop_filter = None #['AFR','EUR','SAS'] ##removed because technically these are not linked to pop but will be filtered anyway and would just open up whole superpop

for pop in pop_filter:
    samples_to_keep = panel.loc[panel['pop'] == pop, 'sample'].tolist()
    ###adding in pop filters
    callset = allel.read_vcf(tmpvcf, fields=['calldata/GT','variants/REF', 'variants/ALT', 'variants/POS'], samples=samples_to_keep) ##select VCF here ##add 'variants/CHROM' in here if needed in the future

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
        'pop_mean_observed_het': obs_het.mean(),
        'expected_het': exp_het,
        'pop_mean_expected_het': exp_het.mean(),
        'difference': diff,
        'pop_mean_diff': diff.mean(),
        'position': allele_pos,
        'ref': ref_alleles,
        'alt': [a[0].decode() if isinstance(a[0], bytes) else str(a[0]) for a in alt_alleles]
    })
    results['population'] = pop
    df_list.append(results)

df_all = pd.concat(df_list, ignore_index=True)

# pd.set_option('display.max_columns', None)
print(df_all.head()) 
reportname = "pop-" + str(pop_filter) + "_report.csv" ##removed super pop
results.to_csv(reportname)