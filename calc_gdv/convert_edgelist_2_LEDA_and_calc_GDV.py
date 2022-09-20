#
#
#
import logging

import click
# #header section
# LEDA.GRAPH
# string
# float
# -2
# #nodes section
# 5
# |{v1}|
# |{v2}|
# |{v3}|
# |{v4}|
# |{v5}|
#
# #edges section
# 5
# 1 2 0 |{}|
# 1 3 0 |{}|
# 2 3 0 |{}|
# 2 4 0 |{}|
# 3 5 0 |{}|


def convert_edgelist_2_LEDA_and_calc_GDV(edgelist,output_dir,count_py_script_path,sep="\t",graph_rep="NAPS",suffix=''):
    print(edgelist)
    p=os.path.join(output_dir, suffix + ".gw")
    output = open(p,'w')
    # edgelist = pd.read_csv(os.path.join(edgelist_dir,f), delimiter="\t",header=None)
    G = nx.readwrite.edgelist.read_edgelist(edgelist, nodetype=str,delimiter=sep)
    if graph_rep != "NAPS":
        G = nx.readwrite.edgelist.read_weighted_edgelist(edgelist, nodetype=str,delimiter=sep)
        G = nx.Graph(((source, target, attr) for source, target, attr in G.edges(data=True) if attr['weight'] > 400))

    output.write("#header section\n")#
    output.write("LEDA.GRAPH\n")
    output.write("string\n")
    output.write("float\n")
    output.write("-2\n")
    output.write("#nodes section\n")
    nodes = list(G.nodes)
    # A5_1A04A
    if graph_rep != "NAPS":
        pass
    else:
        nodes = [n+"_"+suffix for n in nodes]
    print(nodes)
    number_of_nodes = G.number_of_nodes()
    codes_d = dict(zip(nodes,range(1,number_of_nodes+1)))
    print(codes_d)
    output.write(str(number_of_nodes)+"\n")
    for n in nodes:
        output.write(f'|{{{n}}}|\n')
    output.write('\n')
    output.write("#edges section\n")
    number_of_edges = G.number_of_edges()
    output.write(str(number_of_edges)+"\n")
    for n, i in enumerate(G.edges()):
        if graph_rep != "NAPS":
            output.write(f'{codes_d.get(i[0])} {codes_d.get(i[1])} 0 |{{{G.get_edge_data(i[0],i[1])["weight"]}}}|\n')
        else:
            output.write(f'{codes_d.get(i[0]+"_"+suffix)} {codes_d.get(i[1]+"_"+suffix)} 0 |{{}}|\n')
            # print(i)
            # output.write(f'{i[0]} {i[1]} 0 |{{{G.get_edge_data(i[0],i[1])["weight"]}}}|\n')
    output.close()

    # scripts_dir = os.path.join("/groups/pupko/natannag", "consurf_n2v", "huang")
    # script_path =  os.path.join(scripts_dir,'calc_gdv', 'count.py')
    # cmd = f'python {count_py_script_path} {p}'
    cmd = f'cd {CALC_GDV_SCRIPTS}; python ./count.py {p}'
    logging.debug(cmd)
    # subprocess.check_output(cmd,shell=True)
    subprocess.check_output(cmd,shell=True)

def concatenate_GDV_tables(output_dir):
    tables=[]
    for f in os.listdir(output_dir):
        if f.endswith(".ndump2"):
            p= os.path.join(output_dir,f)
            gdv_table = pd.read_csv(p,sep=" ",header=None)
            gdv_table = gdv_table.set_index(0)
            gdv_table = gdv_table.add_prefix("graphlet")
            gdv_table = gdv_table.reset_index()
            gdv_table = gdv_table.rename(columns={0:'pdb_position'})
            tables.append(gdv_table)
    final_GDV_table = pd.concat(tables)
    return final_GDV_table


@click.command()
@click.argument('edgelist', type=click.Path(exists=True))
@click.argument('output_dir', type=click.Path(exists=True))
@click.option('--sep', type=str, default = '\t')
@click.option('--job-title', type=str,default = '')
@click.option('--graph-rep', type=click.Choice(['NAPS', 'PPI'], case_sensitive=True),
              default='NAPS',show_default=True,help='NAPS/PPI if graph represents protein structure or protein interaction network, respectively')
def main(edgelist,output_dir, sep,graph_rep,job_title):
    convert_edgelist_2_LEDA_and_calc_GDV(edgelist,output_dir,COUNT_GDV_SCRIPT,sep,graph_rep,suffix=job_title)
    final_GDV_table= concatenate_GDV_tables(output_dir)
    final_GDV_table.to_csv(os.path.join(output_dir,job_title+'.gdv.csv'),index=False)



if __name__ == '__main__':
    import pandas as pd
    import os
    import subprocess
    import networkx as nx
    from EvoRator.evorator_CONSTANTS import COUNT_GDV_SCRIPT, CALC_GDV_SCRIPTS
    main()
