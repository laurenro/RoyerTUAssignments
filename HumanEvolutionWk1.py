# RoyerTUAssignments
import allel
import pandas as pd
import collections #import Counter
import numpy as np
import matplotlib.pyplot as plt

tmpvcf = 'C:/Users/laure/Documents/BIOL5131/Altai_1000lines.vcf'
# longvcf = 'C:/Users/laure/Documents/BIOL5131/AltaiNea.hg19_1000g.22.mod.vcf.gz' ##update to just 1000lines to run from GitHub repo

# readvcf = allel.read_vcf(tmpvcf)
# keys = readvcf.keys()
# print(keys) 

#############
### update to specific call rather than full df for counts
#############

subset = allel.read_vcf(tmpvcf, fields = ['variants/REF', 'variants/ALT', 'variants/POS'])

alt_alleles = subset['variants/ALT']
ref_alleles = subset['variants/REF']
allele_pos = subset['variants/POS']

###################
combinationcounts = collections.Counter()
varposcounts = collections.Counter()

for ref, alt_list, alpos in zip(ref_alleles, alt_alleles, allele_pos):
    # Flag to check if we found any valid ALT
    has_valid_alt = False
    
    for alt in alt_list:
        if alt and alt != '.':
            combinationcounts[(ref, alt)] +=1
            has_valid_alt = True
            varposcounts[alpos] +=1
    
    # If no valid ALT, count as ref -> ref
    if not has_valid_alt:
        combinationcounts[(ref, ref)] +=1

###Print results
# for (ref, alt), count in combinationcounts.items():
#     print(f"{ref} -> {alt}: {count}") ##turned off so we don't get too much output, turn back on for final!
###Print results in table for just combo
countdf = pd.DataFrame.from_dict(combinationcounts, orient= "index", columns=["count"])
countdf.index =  pd.MultiIndex.from_tuples(countdf.index, names=["REF", "ALT"])
table = countdf["count"].unstack(fill_value=0)

# table.to_csv('tmp_val.csv', index=False) ##turned off so we don't get too much output, turn back on for final!

####print results in table for bp differences
pos_df = pd.DataFrame.from_dict(varposcounts, orient="index", columns=["num_differences"])
pos_df.index.name = "POS"
# pos_df = pos_df.sort_index()
pos_df['binned_pos'] = (pos_df.index // 10) * 10 ###update to 100k for final run
binned = pos_df.groupby("binned_pos")["num_differences"].sum().reset_index()

print("\nDifferences per position:")
print(pos_df)

plt.figure()
plt.bar(binned["binned_pos"], binned["num_differences"])
plt.xlabel("Pos")
plt.ylabel("Diff")
plt.title("Vart Diff/Pos")
plt.show()
plt.savefig('tmp_val.png') ##turned off so we don't get too much output, turn back on for final!