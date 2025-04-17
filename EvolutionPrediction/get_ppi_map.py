import pickle
import sys
import subprocess
import os
import pickle
import numpy as np
import pandas as pd

INPUT = sys.argv[1] # DIR IN WHICH EACH DIR CONTAINS AN EVORATIR OUTPUT DIR
protein_interaction_d = []
i=0
for f in os.listdir(INPUT):
    if i==100:
        break
    tmp_outdir = os.path.join(INPUT,f)
    try:
        with open(os.path.join(tmp_outdir, 'evorator_output_path.txt'), 'r') as f1:
            evorator_df_path = f1.read().strip()
            evorator_df = pd.read_csv(evorator_df_path)
            protein_interaction_d.append(dict(zip(evorator_df.pdb_position, evorator_df.protein_interaction)))
    except Exception as e:
        print(e)
        continue
    i+=1

protein_interaction_d = {pos:ppi for element in protein_interaction_d for pos,ppi in element.items()}
print('saving protein interaction map of size',len(protein_interaction_d),'\n')
print('items',list(protein_interaction_d.items())[:10],'\n')
print('contact','Contact' in list(protein_interaction_d.values()),'\n')
print('interface','Interface' in list(protein_interaction_d.values()),'\n')
with open(os.path.join(INPUT, "ppi_map.data"), 'wb') as f:
    pickle.dump(protein_interaction_d,f)
