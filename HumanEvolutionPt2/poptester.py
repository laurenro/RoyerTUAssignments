###poptester
import allel
import pandas as pd
import numpy as np

tmpvcf = 'C:/Users/laure/RoyerTUAssignments/HumanEvolutionPt2/test_500.vcf.gz'
callset = allel.read_vcf(tmpvcf, fields=["variants/POS", "variants/REF", "variants/ALT", "calldata/GT"])#, samples=True)

###try to read panel
paneltest = pd.read_csv('integrated_call_samples_v3.20130502.ALL.panel', sep='\t')
print(paneltest.head())