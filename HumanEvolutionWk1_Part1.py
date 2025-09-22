import allel
import collections
import pandas as pd

# Load VCF file
vcf_path = 'C:/Users/laure/Documents/BIOL5131/AltaiNea.hg19_1000g.22.mod.vcf.gz'

# Read REF and ALT alleles
callset = allel.read_vcf(vcf_path, fields=['variants/REF', 'variants/ALT'])

refs = callset['variants/REF']
alts = callset['variants/ALT']

# Count combinations
counts = collections.Counter()

for ref, alt_list in zip(refs, alts):
    has_valid_alt = False
    for alt in alt_list:
        if alt and alt != '.':
            counts[(ref, alt)] += 1
            has_valid_alt = True
    if not has_valid_alt:
        counts[(ref, ref)] += 1

# Build DataFrame (matrix form)
df = pd.DataFrame.from_dict(counts, orient="index", columns=["count"])
df.index = pd.MultiIndex.from_tuples(df.index, names=["REF", "ALT"])
table = df["count"].unstack(fill_value=0)

table.to_csv('cleaned-up_fina_sub_part1.csv')