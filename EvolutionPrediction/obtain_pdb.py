import networkx as nx
import sys
import numpy as np
import os
import matplotlib.pyplot as plt
import click

laptop = False if (os.path.exists("/groups/pupko/natannag") or os.path.exists("/bioseq")) else True

if laptop:
    path2github = "../"
else:
    path2github = "/groups/pupko/natannag"

sys.path.append(path2github)

if laptop:
    path2proteinnet = "C:/Users/natan/Documents/ScanNet_dev"
else:
    path2proteinnet = '/groups/pupko/natannag/ScanNet_dev/'
sys.path.append(path2proteinnet)

from preprocessing import PDBio,PDB_processing
def obtain_pdb(identifier,job_dir,all_chains_flg = False):
    ## identifier = 1a3x_0-A / 1a3x / list of such?
    ## job_dir = where structures are saved
    ## all_chains_flg = set to True/False to generate a graph for the biological unit/a single chain extracted from the biological unit

    assert (len(identifier)==4 and  all_chains_flg!=False) or (len(identifier)>4 and  all_chains_flg!=True) , f"pdb identifier {identifier} format is not compatible with all_chains_flg = {all_chains_flg}. set all_chains_flg to True or specify the pdb_id "
    file, chain_ids = PDBio.getPDB(identifier, biounit=True,structures_folder=job_dir)
    final_pdb_file = os.path.join(''.join(os.path.split(file)[:-1]), os.path.split(file)[-1].split(".")[0] + '.pdb')
    chains = PDBio.load_chains(file=file, chain_ids=chain_ids)[1]

    sequence_from_pdb = PDB_processing.process_chain(chains)[0]
    backbone_coordinates = PDB_processing.process_chain(chains)[2]
    residue_indices = PDB_processing.get_PDB_indices(chains, return_model=True, return_chain=True)

    if all_chains_flg:
        print(chain_ids)
        if chain_ids == 'all':
            chain_ids = [(0, i.id) for i in chains]
        PDBio.extract_chains(file, [*chain_ids], final_pdb_file)
    else:
        final_pdb_file = file


    return final_pdb_file, chain_ids, sequence_from_pdb, residue_indices, backbone_coordinates
