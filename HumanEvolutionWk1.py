# RoyerTUAssignments
import allel
import pandas as pd
import collections
import matplotlib.pyplot as plt

# tmpvcf = 'Altai_1000lines.vcf'
longvcf = 'C:/Users/laure/Documents/BIOL5131/AltaiNea.hg19_1000g.22.mod.vcf.gz' ##update to just 1000lines to run from GitHub repo

# # readvcf = allel.read_vcf(longvcf)
# # keys = readvcf.keys()
# # print(keys) 

subset = allel.read_vcf(longvcf, fields = ['variants/REF', 'variants/ALT', 'variants/POS'])

alt_alleles = subset['variants/ALT']
ref_alleles = subset['variants/REF']
allele_pos = subset['variants/POS']

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

###Print results in table for just combo
countdf = pd.DataFrame.from_dict(combinationcounts, orient= "index", columns=["count"])
countdf.index =  pd.MultiIndex.from_tuples(countdf.index, names=["REF", "ALT"])
table = countdf["count"].unstack(fill_value=0)
alleles = ['A', 'C', 'G', 'T']
table = table.reindex(index=alleles,columns=alleles,fill_value=0)
table.to_csv('simple-allele-count.csv') ################ part 1 table

####print results in table for bp differences
pos_df = pd.DataFrame.from_dict(varposcounts, orient="index", columns=["num_differences"])
pos_df.index.name = "POS"
# pos_df = pos_df.sort_index() ######################## part 2 table
pos_df['binned_pos'] = (pos_df.index // 100000) * 100000
binned = pos_df.groupby("binned_pos")["num_differences"].sum().reset_index()
binned.to_csv('long_pos.csv', index = False) ##turned off so we don't get too much output, turn back on for final!

# print("\nDifferences per position by 100k:")
# print(pos_df)

####################### part 3 plot and save
plt.figure()
plt.bar(binned["binned_pos"], binned["num_differences"], width=80000, align='center')
plt.xlabel("Pos")
plt.ylabel("Diff")
plt.title("Vart Diff/Pos")
# plt.show()
plt.ylim(0, max(2, binned["num_differences"].max() * 1.1))
plt.tight_layout()
plt.savefig('long_val.png', dpi=300)  ##turned off so we don't get too much output, turn back on for final!
plt.close()