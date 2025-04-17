import numpy as np
from Bio.PDB import Selection
from itertools import chain

import sys

from Bio.PDB import PDBParser
import os
import pickle

laptop = False if os.path.exists("/groups/pupko/natannag/natan_git") else True
if laptop:
    path2github = "../"
    path2scannet = "C:/Users/natan/Documents/ScanNet_dev/"
    # path2scannet = "C:/Users/natan/Documents/ScanNet/"
    # path2scannet = "C:/Users/natan/Documents/missense_pathogenecity/ScanNet-main/ScanNet/"""
    path2evolutionprediction = "C:/Users/natan/Documents/EvolutionPrediction/"
else:
    path2scannet = "/groups/pupko/natannag/natan_git/ScanNet_dev"
    path2github = "../"


def is_residue(residue):
    try:
        return (residue.get_id()[0] in hetresidue_field) & (residue.resname in residue_dictionary.keys())
    except:
        return False


residue_dictionary = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
                      'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
                      'GLY': 'G', 'HIS': 'H', 'HSD': 'H', 'HSE': 'H',
                      'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
                      'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M',
                      'MSE': 'M',
                      'PTR': 'Y',
                      'TYS': 'Y',
                      'SEP': 'S',
                      'TPO': 'T',
                      'HIP': 'H', }

hetresidue_field = [' '] + ['H_%s' % name for name in residue_dictionary.keys()]

sys.path.append(path2github)
proteingym_PDB = 'proteingym_PDB/pdb_files'

residue_pdb_index = []
mean_pLDDT_map = {}
for f in os.listdir(proteingym_PDB):
    parser = PDBParser(QUIET=True)
    if f.endswith(".pdb"):
        print(f)
        with open(os.path.join(proteingym_PDB, f), 'r') as pdb_file:
            structure = parser.get_structure('struct', pdb_file)
            chain_list = [c for c in structure.get_chains()]
            sequence = ''
            bfactors = []
            for residue in Selection.unfold_entities(chain_list[0],
                                                     'R'):  # Part 2: Data Acquisition - In section 2.3, it is uncear from which chain to extract the sequence
                if is_residue(residue):
                    sequence += residue_dictionary[residue.resname]
                    residue_pdb_index.append(residue.id[1])
                    for atom in Selection.unfold_entities(residue, 'A'):
                        if atom.get_id() == 'CA':
                            # print(atom.get_bfactor())
                            bfactors.append(atom.get_bfactor())
                            break
            # print(len(sequence))
            # print(len(bfactors))
            # print(sum(bfactors))
            mean_pLDDT = (sum(bfactors)) / len(sequence)  # print(Selection.get_unique_parents([atom]))

        mean_pLDDT_map[f.split(".")[0].split("_")[-1]] = mean_pLDDT

with open('mean_pLDDT_map.data', 'wb') as mean_pLDDT_map_file:
    pickle.dump(mean_pLDDT_map, mean_pLDDT_map_file)
