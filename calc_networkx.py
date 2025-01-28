import click

def get_list_of_graphs(edgelists_dir,sep=None):
    list_of_graphs = []
    list_of_node_names_map = []
    for j, f in enumerate(os.listdir(edgelists_dir)):
        try:

            print(j, f)
            suffix = f.split("_")[0]

            edgelist = os.path.join(edgelists_dir, f)
            G = nx.readwrite.edgelist.read_weighted_edgelist(edgelist, nodetype=str,delimiter=sep)
            list_of_graphs.append((G,suffix))
        except:
            continue

    return list_of_graphs

def rename_nodes(list_of_graphs):
    list_of_node_names_map = [dict(zip(i[0].nodes, [n + "_" + i[1] for n in list(i[0].nodes)])) for i in list_of_graphs]
    names_map = get_nodes_d(list_of_node_names_map)
    LG = [nx.relabel_nodes(i[0], names_map) for i in list_of_graphs]
    return LG

def get_nodes_d(list_of_dicts):
    super_dict = {}
    for d in list_of_dicts:
        for k, v in d.items():  # d.items() in Python 3+
            super_dict[k]=v
    return super_dict

def get_neighbors_4_each_node(G,suffix=""):

    neighbors_dict = {str(i)+suffix:[str(n)+suffix for n in G.neighbors(i)] for i in list(G.nodes)}

    return neighbors_dict

#Centrality
def get_eigenvector_centrality(G,suffix =''):
    print(G)
    nodes_eigenvector_centrality_dicts = []
    feature_d = {}
    try:
        centrality = nx.eigenvector_centrality(G,max_iter=1000)
    except:
        centrality = nx.eigenvector_centrality(G,max_iter=10000)
    centrality = {k + suffix : v for k,v in centrality.items()}
    feature_d['eigenvector_centrality'] = centrality
    final = pd.DataFrame.from_dict(feature_d, orient='columns')
    return final

#Centrality
def get_betweenness_centrality(G,suffix=""):
    list_of_dicts = []
    feature_d = {}

    centrality = nx.betweenness_centrality(G)
    centrality = {k + suffix : v for k,v in centrality.items()}
    list_of_dicts.append(centrality)

    nodes_d=get_nodes_d(list_of_dicts)
    feature_d['betweenness_centrality'] = nodes_d
    final = pd.DataFrame.from_dict(feature_d, orient='columns')
    return final

def get_closeness_centrality(G):
    list_of_dicts = []
    feature_d = {}

    for G,suffix in G:
        centrality = nx.closeness_centrality(G)
        centrality = {k + "_" + suffix : v for k,v in centrality.items()}
        list_of_dicts.append(centrality)

    nodes_d=get_nodes_d(list_of_dicts)
    feature_d['closeness_centrality'] = nodes_d
    final = pd.DataFrame.from_dict(feature_d, orient='columns')
    return final



#degree
def get_node_degree(G,suffix=""):
    list_of_dicts = []
    centrality = dict(nx.degree(G))
    centrality = {k + suffix : v for k,v in centrality.items()}
    list_of_dicts.append(centrality)
    feature_d = {}
    nodes_d=get_nodes_d(list_of_dicts)
    feature_d['node_degree'] = nodes_d
    final = pd.DataFrame.from_dict(feature_d, orient='columns')
    return final

# Centrality
def degree_centrality(G,suffix=""):
    list_of_dicts = []
    centrality = dict(nx.degree_centrality(G))
    print(centrality)
    centrality = {k + suffix : v for k,v in centrality.items()}
    list_of_dicts.append(centrality)
    feature_d = {}
    nodes_d=get_nodes_d(list_of_dicts)
    feature_d['degree_centrality'] = nodes_d
    final = pd.DataFrame.from_dict(feature_d, orient='columns')
    return final

# Assortativity
def get_average_neighbor_degree(G,suffix=""):
    list_of_dicts = []
    average_neighbor_degree = nx.average_neighbor_degree(G,weight='weight')
    print(average_neighbor_degree)
    average_neighbor_degree = {k + suffix : v for k,v in average_neighbor_degree.items()}
    list_of_dicts.append(average_neighbor_degree)
    feature_d = {}
    nodes_d=get_nodes_d(list_of_dicts)
    feature_d['average_neighbor_degree'] = nodes_d
    final = pd.DataFrame.from_dict(feature_d, orient='columns')
    return final

# Clustering
def get_clustering_coefficient(G,suffix=''):
    list_of_dicts = []
    clustering_c = nx.clustering(G,weight='weight')
    print(clustering_c)
    clustering_c = {k + suffix : v for k,v in clustering_c.items()}
    list_of_dicts.append(clustering_c)
    feature_d = {}
    nodes_d=get_nodes_d(list_of_dicts)
    feature_d['clustering_coefficient'] = nodes_d
    final = pd.DataFrame.from_dict(feature_d, orient='columns')
    return final

# Clique
def get_k_cliques(G,suffix=""):
    clique_dicts = []
    clique_d = {}
    cliques = nx.find_cliques(G)
    c = list(cliques)
    for i in c:
        for n,j in enumerate(i):
            i[n] = i[n] + suffix

    for i in c:
        for j in i:
            if j in clique_d:
                    clique_d[j].append(str(len(i))+"_clique")
            else:
                    clique_d[j] = [str(len(i))+"_clique"]

    for k in clique_d:
        clique_d[k] = dict(Counter(clique_d[k]))

    clique_dicts.append(clique_d)


    nodes_d=get_nodes_d(clique_dicts)
    feature_d={}
    feature_d['k_cliques'] = nodes_d
    final = pd.DataFrame.from_dict(nodes_d, orient='index').fillna(0)
    print(final)
    return final

@click.command()
@click.argument('edgelist', type=click.Path(exists=True))
@click.argument('calc_networkx_output_dir', type=click.Path(exists=True))
@click.option('--job-title',type=str,default='',show_default=True,help='Insert job title')
def main(edgelist, calc_networkx_output_dir,job_title):

    out = os.path.join(calc_networkx_output_dir,f'{job_title}.naps.csv')
    G = nx.readwrite.edgelist.read_weighted_edgelist(edgelist, nodetype=str, delimiter='\t')
    neighbors_dict = get_neighbors_4_each_node(G,suffix="_"+job_title)
    d1=get_eigenvector_centrality(G,suffix="_"+job_title)
    d2=get_betweenness_centrality(G,suffix="_"+job_title)
    d3=get_node_degree(G,suffix="_"+job_title)
    d4=degree_centrality(G,suffix="_"+job_title)
    d5=get_average_neighbor_degree(G,suffix="_"+job_title)
    d6=get_clustering_coefficient(G,suffix="_"+job_title)
    d7=get_k_cliques(G,suffix="_"+job_title)
    neighbors_df = pd.DataFrame.from_dict(neighbors_dict, orient='index').add_prefix("neighbor_pos_")
    final_table = pd.concat([d1,d2,d3,d4,d5,d6,d7,neighbors_df],axis=1)
    final_table = final_table.reset_index()
    final_table = final_table.rename(columns={'index':'pdb_position'})
    final_table.to_csv(out,index=False)



if __name__ == '__main__':
    import time
    import numpy as np
    import pandas as pd
    import subprocess
    import re
    import networkx as nx
    from collections import Counter
    import os
    from networkx.algorithms.community import k_clique_communities
    main()
