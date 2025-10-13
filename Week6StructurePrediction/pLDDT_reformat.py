import os, sys, glob
import json
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

VERBOSE = True


### Load in the predicted pLDDT values from AF2

# AF2_dir = 'C:/Users/laure/Downloads/ChimeraX/AlphaFold/prediction_1'
# json_files = glob.glob( os.path.join(AF2_dir,'af168_scores_rank_00?_alphafold2_ptm_model_?_seed_000.json'))
# json_files.sort()   # sort in place
# n = len(json_files)
# af2_labels = [f'AF2 model{i}' for i in range(n)]
# print('json_files', json_files)

AF3_dir = 'C:/Users/laure/Downloads/fold_2025_10_12_16_58'
csv_files = glob.glob( os.path.join(AF3_dir,'model?_test_per_residue_plddt.csv'))
csv_files.sort()   # sort in place
n = len(csv_files)
af3_labels = [f'AF3 model{i}' for i in range(n)]
print('json_files', csv_files)

###this is for both
all_data  = []

# # Open and read the JSON file ###THIS IS FOR AF2
# for json_file in json_files:
#     with open(json_file, 'r') as file:
#         data = json.load(file)
#     if VERBOSE:
#         print(data["plddt"])
#     all_data.append( data["plddt"] )   # these are per-residue c

# Open and read the internally created CSV file w per resiude info ###THIS IS FOR AF3
for csv_file in csv_files:
    with open(csv_file, 'r') as file:
        data = pd.read_csv(file)
    if VERBOSE:
        print(data["plddt"])
    all_data.append( data["plddt"] )   # these are per-residue c

all_data = np.array(all_data).transpose()
print('all_data.shape', all_data.shape)

### Alternatively, we could read in the pLDDT values from the `best_model.pdb' B-factors
if (0):
    # NOTE: these pLDDT values are per-atom, not per-residue, so we need to have a mapping from
    #       atomindex to residue index
    import subprocess
    # cmd = f'cat {os.path.join(AF2_dir,"best_model.pdb")} | grep ATOM' ###UPDATE DIR HERE FOR AF2/3 PLOTS
    cmd = f'cat {os.path.join(AF3_dir,"fold_2025_10_12_16_58_model_0.pdb")} | grep ATOM' ###UPDATE DIR HERE FOR AF2/3 PLOTS, THIS AF3 ONLY WORKS AFTER MANUAL PDB CONV
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    print('result', result)
    
    lines = result.stdout.split('\n')
    chain_of_atomnum, resnum_of_atomnum = {}, {}
    atomnums, chains, resnums = [], [], []
    for line in lines:
        try:
            fields = line.split()   # ATOM   3244  CD2 HIS B 209   
            atomnum, chain, resnum = int(fields[1]), fields[4], int(fields[5])
            chain_of_atomnum[atomnum] = chain
            resnum_of_atomnum[atomnum] = resnum
            atomnums.append(atomnum)
            chains.append(chain)
            resnums.append(resnum)
        except:
            pass

### Load the ACTUAL lddt from a csv  (obtained from https://swissmodel.expasy.org/assess )
actual_lddt_csv = 'C:/Users/laure/Downloads/AF3_model1_3qVAok_01_lddt.csv' ###UPDATE CSV HERE FOR AF 2/3
actual = pd.read_csv(actual_lddt_csv)
print(actual.columns)
print(actual["lddt"])
print(actual["mdl_res_no"])


plt.figure(figsize=(10.,4.))
plt.plot(actual["mdl_res_no"], all_data/100.0, '.', label=af3_labels, ms=1)    # divide by 100 to convert from percent ####CHANGE AF2/AF3 HERE
plt.plot(actual['mdl_res_no'], actual["lddt"], label='actual LDDT')
plt.xlabel('model residue number')
plt.ylabel('pLDDT / LDDT')
plt.legend(loc='best')
plt.tight_layout()
plt.show()

# outpng = 'AF3_Model1_pLDDT_vs_actual.png' ###UPDATE NAME HERE FOR AF 2/3
# plt.savefig(outpng)
# print('Wrote:', outpng)


