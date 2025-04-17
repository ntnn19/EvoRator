from Bio.PDB import PDBParser
import os
# model_url = f'https://alphafold.ebi.ac.uk/files/{alphafold_ID}-model_{database_version}.pdb'
import json
import pickle
output = 'uniprot_to_pae_d.pkl'
d={}
model_path = r'C:\Users\natan\PDB'
for file in os.listdir(model_path):
    if file.startswith('AF'):
        uniprot_id = file.split("_")[-1].split(".")[0]
        model_pae_url = f'https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-predicted_aligned_error_v4.json'
        os.system(f'curl {model_pae_url} -o proteingym_PDB/pae_files/{uniprot_id}_pae.json')
        j = json.load(open(f'proteingym_PDB/pae_files/{uniprot_id}_pae.json','r'))
        d[uniprot_id]=np.asarray(j[0]['predicted_aligned_error']).flatten().mean()

print(d)
with open(output,'wb') as f:
    pickle.dump(d,f)