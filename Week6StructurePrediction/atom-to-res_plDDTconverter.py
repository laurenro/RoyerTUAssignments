import json
import numpy as np
import csv
import matplotlib.pyplot as plt
from Bio.PDB import PDBParser

# === User Inputs ===
pdb_file = 'C:/Users/laure/Downloads/fold_2025_10_12_16_58/fold_2025_10_12_16_58_model_4.pdb'                 # AlphaFold3 predicted structure
json_file = 'C:/Users/laure/Downloads/fold_2025_10_12_16_58/fold_2025_10_12_16_58_full_data_4.json'   # AlphaFold3 output JSON with atom_plddts
output_csv = 'model4_test_per_residue_plddt.csv'

# === Load atom_plddts from AlphaFold3 JSON ===
with open(json_file, 'r') as f:
    af_json = json.load(f)
atom_plddts = af_json['atom_plddts']

# === Parse the PDB structure ===
parser = PDBParser(QUIET=True)
structure = parser.get_structure("AF3", pdb_file)

residue_plddt_dict = {}  # (chain_id, resnum) -> list of atom plddts
residue_info = {}        # (chain_id, resnum) -> residue_name

atom_idx = 0

for model in structure:
    for chain in model:
        for residue in chain:
            res_id = residue.id
            if res_id[0] != " ":  # Skip hetero/water
                continue
            chain_id = chain.id
            res_num = res_id[1]
            res_key = (chain_id, res_num)

            atom_plddts_list = []
            for atom in residue:
                if atom_idx >= len(atom_plddts):
                    raise ValueError("Mismatch: More atoms in PDB than pLDDTs in JSON.")
                atom_plddt = atom_plddts[atom_idx]
                atom_plddts_list.append(atom_plddt)
                atom_idx += 1

            residue_plddt_dict[res_key] = atom_plddts_list
            residue_info[res_key] = residue.get_resname()

# === Compute per-residue pLDDT (mean of atom values) ===
results = []
for (chain_id, resnum), atom_scores in sorted(residue_plddt_dict.items(), key=lambda x: x[0][1]):
    mean_plddt = np.mean(atom_scores)
    res_name = residue_info[(chain_id, resnum)]
    results.append({
        'chain_id': chain_id,
        'residue_number': resnum,
        'residue_name': res_name,
        'plddt': mean_plddt
    })

# === Save to CSV ===
with open(output_csv, 'w', newline='') as csvfile:
    fieldnames = ['chain_id', 'residue_number', 'residue_name', 'plddt']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()
    for row in results:
        writer.writerow(row)

print(f"✅ Per-residue pLDDT saved to: {output_csv}")

# === Optional: Plot pLDDT per residue ===
residue_numbers = [row['residue_number'] for row in results]
plddts = [row['plddt'] for row in results]

plt.figure(figsize=(12, 5))
plt.plot(residue_numbers, plddts, label='AlphaFold3 per-residue pLDDT', color='blue')
plt.xlabel('Residue Number')
plt.ylabel('pLDDT Score')
plt.title('AlphaFold3: Per-residue Predicted LDDT (pLDDT)')
plt.ylim(0, 100)
plt.grid(True)
plt.tight_layout()
plt.legend()
plt.show()
