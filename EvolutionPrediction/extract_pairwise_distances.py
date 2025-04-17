import networkx as nx
import sys
import numpy as np
import os
import obtain_pdb
import matplotlib.pyplot as plt
import click

laptop = False if os.path.exists("/groups/pupko/natannag") else True

if laptop:
    path2github = "../"
else:
    path2github = "/groups/pupko/natannag"

sys.path.append(path2github)

def get_edge_list(identifier,backbone_coordinates,residue_pdb_index,job_dir):
    output_name_all = identifier.split("_")[0].upper()+"_edgelist.txt"
    output_all = os.path.join(job_dir,output_name_all)
    Calpha_coordinates = backbone_coordinates[:, 2, :]
    Calpha_Calpha_distance = np.sqrt(
        ((Calpha_coordinates[np.newaxis] - Calpha_coordinates[:, np.newaxis]) ** 2).sum(-1))
    lower_cutoff = 0
    upper_cutoff = 7
    Calpha_Calpha_distance_filtered = (float(lower_cutoff) <= Calpha_Calpha_distance) & (
                Calpha_Calpha_distance <= float(upper_cutoff))
    Calpha_Calpha_distance_filtered = Calpha_Calpha_distance_filtered.astype(int)
    print(Calpha_Calpha_distance_filtered.shape)
    plt.savefig(os.path.join(job_dir,'debug_naps.png'))
    residue_pdb_index_with_chain_id = []

    for i in residue_pdb_index:
        residue_pdb_index_with_chain_id.append(str(i[1])+str(i[-1]))
    int2label = dict(zip(list(range(Calpha_Calpha_distance_filtered.shape[0])), residue_pdb_index_with_chain_id))
    G = nx.from_numpy_matrix(Calpha_Calpha_distance_filtered)
    G = nx.relabel_nodes(G, int2label)
    G.remove_edges_from(nx.selfloop_edges(G))
    print(G.nodes)


    nx.write_edgelist(G,output_all , data=False,delimiter='\t')
    return output_all



