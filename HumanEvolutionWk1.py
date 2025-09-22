# RoyerTUAssignments
import allel
import pandas as pd
import collections #import Counter
import numpy as np

tmpvcf = 'C:/Users/laure/Documents/BIOL5131/Altai_1000lines.vcf'
# longvcf = 'C:/Users/laure/Documents/BIOL5131/AltaiNea.hg19_1000g.22.mod.vcf.gz'

# readvcf = allel.read_vcf(tmpvcf)
# keys = readvcf.keys()
# print(keys)

#############
### update to specific call rather than full df for counts
#############

subset = allel.read_vcf(tmpvcf, fields = ['variants/REF', 'variants/ALT'])

alt_alleles = subset['variants/ALT']
ref_alleles = subset['variants/REF']

###################
combinationcounts = collections.Counter()

for ref, alt_list in zip(ref_alleles, alt_alleles):
    # Flag to check if we found any valid ALT
    has_valid_alt = False
    
    for alt in alt_list:
        if alt and alt != '.':
            combinationcounts[(ref, alt)] += 1
            has_valid_alt = True
    
    # If no valid ALT, count as ref -> ref
    if not has_valid_alt:
        combinationcounts[(ref, ref)] += 1

# Print results
for (ref, alt), count in combinationcounts.items():
    print(f"{ref} -> {alt}: {count}")