import evorator_CONSTANTS as CONSTS  # from /bioseq/natan_conservation_webserver/auxiliaries/
import string
import itertools
import time
import pandas as pd
import numpy as np
import os
import subprocess
from Bio.PDB import PDBParser
from functools import reduce
import logging
import sys



#
# import click
# @click.command()
# @click.argument('edgelist', type=click.Path(exists=True))
# @click.argument('pdb_file', type=click.Path(exists=True))
# @click.argument('pdb_chain', type=str)
# # @click.argument('calc_wcn_output_dir', type=click.Path(exists=False))
# # @click.argument('dssp_output_dir', type=click.Path(exists=True))
# # @click.argument('calc_rsa_output_dir',type=click.Path(exists=False))
# # @click.argument('calc_networkx_output',type=click.Path(exists=False))
# # @click.argument('calc_embedding_output',type=click.Path(exists=False))
# # @click.argument('calc_gdv_output',type=click.Path(exists=False))
# @click.argument('catalytic_sites',type=click.Path(exists=True))
# # @click.argument('disordered_sites',type=click.Path(exists=True))
# # @click.argument('binding_sites',type=click.Path(exists=True))
# # @click.argument('rates_table_file',type=click.Path(exists=True))
# @click.argument('results_dir', type=click.Path(exists=True))
# @click.option('--job-title',type=str,default='',show_default=True,help='Insert job title')
def get_neigh_features(prop_d, curr_prop, df_merged, neigh_d):
    for i, r in enumerate(curr_prop):
        if 'total_' + str(r) + '_neigh' in df_merged.columns.tolist() or str(r)=='B':
            col_name = 'total_ss_' + str(r) + '_neigh'
        else:
            col_name = 'total_' + str(r) + '_neigh'
        df_merged[col_name] = [[prop_d[k] for k in neigh_d[row] if prop_d.get(k, 0) == r] for
                                                   row in neigh_d if row in df_merged.pdb_position.tolist()]
        df_merged[col_name] = df_merged[col_name].apply(len)

def position2property_dict(df_merged,var):
    try:
        return dict(zip(df_merged['pdb_position'], df_merged[var]))
    except:
        return dict(zip(df_merged['pdb_position'], [np.nan]*df_merged.shape[0]))

def extract_features(edgelist, pdb_file_complete, pdb_file_single_chain, pdb_chain, catalytic_sites, results_dir,
                     consurf_output='', job_title=''):
    ''' This script extract features from a given PDB_FILE , PDB_CHAIN pair and an unweighted EDGELIST file.
    Returns n x m matrix where n is the number of positions in the PDB_FILE and m is the number of features'''

    job_result_dir = results_dir
    if os.path.exists("/groups/pupko/natannag/"):
        scripts_dir = CONSTS.EVORATOR_EXEC
        sys.path.append(os.path.join(scripts_dir,"glycosylator-master"))
        try:
            import glycosylator as gl
        except:
            pass

    else:
        scripts_dir = ""

    #    if job_title== '':
    #       job_title = os.path.split(pdb_file_single_chain)[-1].split(".")[0]


    # run DSSP
    dssp_output = os.path.join(job_result_dir, job_title + ".dssp")
    try:
        if not os.path.exists(dssp_output):
            logging.debug('Running DSSP')
            # cmd = f'module unload python/python-anaconda3.7-itaym;module load dssp-4.0; module load boost/boost-1-75-0; module load gcc/gcc-8.2.0; mkdssp -i {pdb_file_single_chain} --output-format=dssp -o {dssp_output}'

            with open(pdb_file_complete,'r') as f:
                content = f.read()
            if 'CRYST1' not in content:
                with open(pdb_file_complete, 'w') as f:
                    f.write('CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1           1          \n')
                    f.write(content)
            cmd = f'module unload python/python-anaconda3.7-itaym;module load dssp-4.0; module load boost/boost-1-75-0; module load gcc/gcc-8.2.0; mkdssp -i {pdb_file_complete} --output-format=dssp -o {dssp_output}'
            logging.debug(cmd)
            subprocess.check_output(cmd, shell=True)
    except:
        logging.debug("DSSP failed")
        pass

    # subprocess.check_output(f'module load python/python-anaconda3.7-itaym', shell=True)

    # calc RSA
    try:
        if not os.path.exists(dssp_output + ".rsa.csv"):
            logging.debug('Running calc_rsa')
            calc_rsa_script = os.path.join(scripts_dir, "calc_rsa.py")
            cmd = f'python {calc_rsa_script} {dssp_output} {job_result_dir}'
            logging.debug(cmd)
            subprocess.check_output(cmd, shell=True)

    except:
        logging.debug("DSSP output file was not found")
        pass
        # raise FileNotFoundError

    # calc WCN
    try:
        if not os.path.exists(
                os.path.join(job_result_dir, os.path.split(pdb_file_complete)[-1].split(".")[0] + ".wcn.csv")):
            logging.debug('Running calc_wcn')
            calc_wcn_script = os.path.join(scripts_dir, "calc_wcn.py")
            cmd = f'python {calc_wcn_script} {pdb_file_complete} -d {job_result_dir}'
            logging.debug(cmd)
            subprocess.check_output(cmd, shell=True)
    except:
        logging.debug("WCN cannot be calculated")
        pass
    # calc NETWORKX FEATURES
    try:
        if not os.path.exists(os.path.join(job_result_dir, job_title + ".naps.csv")):
            logging.debug('Running calc_networkx')
            calc_networkx_script = os.path.join(scripts_dir, "calc_networkx.py")
            cmd = f'python {calc_networkx_script} {edgelist} {job_result_dir} --job-title={job_title}'
            logging.debug(cmd)
            subprocess.check_output(cmd, shell=True)
    except:
        logging.debug("Edgelist file was not found")
        raise FileNotFoundError

    # calc GDV
    if not os.path.exists(os.path.join(job_result_dir, job_title + ".gdv.csv")):
        logging.debug('Running convert_edgelist_2_LEDA_and_calc_GDV')


        cmd = f'python {CONSTS.CONVERT2LEDA_GDV_SCRIPT} {edgelist} {job_result_dir} --job-title={job_title}'
        logging.debug(cmd)
        subprocess.check_output(cmd, shell=True)

    # get disordered sites, binding sites, and catalytic sites
    if not os.path.exists(os.path.join(job_result_dir, job_title + ".disorder.txt")) or not os.path.exists(
            os.path.join(job_result_dir, job_title + ".binding.txt")):
        logging.debug('Running get_disordered_and_binding_sites.py')
        get_disordered_and_binding_sites_script = os.path.join(scripts_dir, 'get_disordered_and_binding_sites.py')

        # cmd = f'python {get_disordered_and_binding_sites_script} {pdb_file_single_chain} {pdb_chain} {job_result_dir} --job-title={job_title}'
        cmd = f'python {get_disordered_and_binding_sites_script} {pdb_file_complete} {pdb_chain} {job_result_dir} --job-title={job_title}'
        logging.debug(cmd)
        subprocess.check_output(cmd, shell=True)

    # calc node2vec
    # if consurf_output:
    #     if not os.path.exists(os.path.join(job_result_dir, job_title + ".n2v.csv")):
    #         logging.debug('Running run_n2v_cluster.py')
    #         calc_node2vec_script = os.path.join(scripts_dir, 'run_n2v_cluster.py')
    #
    #         cmd = f'python {calc_node2vec_script} {edgelist} {job_result_dir} --job-title={job_title}'
    #         logging.debug(cmd)
    #         subprocess.check_output(cmd, shell=True)

    # get contact residues
    GetInterfaces_root_outdir = os.path.join(job_result_dir, 'GetInterfaces_out')
    if not os.path.exists(GetInterfaces_root_outdir):
        os.makedirs(GetInterfaces_root_outdir,exist_ok=True)
    pdb_code = os.path.split(pdb_file_complete)[-1].split(".")[0]
    logging.debug(f'pdb_code1: {pdb_code}')
    pdb_parser = PDBParser().get_structure(pdb_code, pdb_file_complete)
    logging.debug(f'pdb_code2: {pdb_code}')
    chain_list = [c.get_id() for c in pdb_parser.get_chains()]
    interfaces_pairs = list(itertools.combinations(chain_list, 2))
    # only keep pairs that contain the selected chain
    interfaces_pairs = [i for i in interfaces_pairs if i[0] == pdb_chain or i[1] == pdb_chain]
    logging.debug(f'chain pairs: {interfaces_pairs}')
    if len(chain_list) > 1:
        logging.debug('Running GetInterfaces.py')
        GetInterfaces_script = os.path.join(scripts_dir, 'GetInterfaces.py')
        for i in interfaces_pairs:
            c1 = i[0]
            c2 = i[1]
            if not os.path.exists(os.path.join(GetInterfaces_root_outdir, "out_" + job_title + "_" + c1 + "_" + c2)):
                cmd = f'python {GetInterfaces_script} --f1 {pdb_file_complete} --f2 {pdb_file_complete} --c1 {c1} --c2 {c2} --c 4.5 --i 10 --jobid {job_title + "_" + c1 + "_" + c2} --jobdir {GetInterfaces_root_outdir}'
                logging.debug(cmd)
                try:
                    subprocess.check_output(cmd, shell=True)
                except:
                    continue

        dfs = []
        query = f'molecule_{pdb_chain}.txt'
        for root, subdirectories, files in os.walk(GetInterfaces_root_outdir):
            for file in files:
                subject = os.path.join(root, file)
                if query in subject and os.path.getsize(subject) != 0:
                    contacting_residues_table = pd.read_csv(subject, header=None, sep=' ')
                    dfs.append(contacting_residues_table)
        try:
            contacting_residues_table = pd.concat(dfs, axis=0)
        # A468_2PZR

            contacting_residues_table['pdb_position'] = contacting_residues_table[0] + contacting_residues_table.astype(str)[1] + "_" + job_title
            logging.debug(contacting_residues_table)
        except:
            contacting_residues_table = pd.DataFrame(np.nan, index=list(range(2000)), columns=['pdb_position',1,2])
            logging.debug(contacting_residues_table)

    # get glycosylated residues using glycosylator
    myGlycosylator = gl.Glycosylator(os.path.join(gl.GLYCOSYLATOR_PATH, 'support/toppar_charmm/carbohydrates.rtf'),
                                     os.path.join(gl.GLYCOSYLATOR_PATH, 'support/toppar_charmm/carbohydrates.prm'))

    try:
        myGlycosylator.load_glycoprotein(pdb_file_complete)
    # myGlycosylator.load_glycoprotein(pdb_file_single_chain)
        glycosylated_positions = list(myGlycosylator.glycans.keys())
        glycosylated_positions = ["".join(k.split(',')) + "_" + job_title for k in glycosylated_positions]
        logging.debug(f'glycosylated: {glycosylated_positions}')
    except:
        glycosylated_positions = []
        logging.debug(f'glycosylated: {glycosylated_positions}')


    dfs = []
    for f in os.listdir(job_result_dir):
        try:
            if f.split(".")[-2] in ['gdv', 'n2v', 'naps', 'rsa', 'wcn'] and 'temp' not in f and '.csv' in f:
                p = os.path.join(job_result_dir, f)

                logging.debug(p)
                df = pd.read_csv(p)
                if f.split(".")[-2] in ['rsa', 'wcn']:
                    df['pdb_position'] = df.chain.astype(str) + df.pdb_position.astype(str) + "_" + job_title
                    df = df[df.chain==pdb_chain]
                dfs.append(df)
        except:
            continue
            # elif f.split(".")[-1] in ['ER','OU','cmds','html','pbs','log','ndump2','pdb','gw','dssp'] or f.split(".")[-2] in ['disorder','binding','cgi_debug'] or f.split("_")[-1] in ['edgelist.txt'] :
            #   continue
            # else:
            #   p = os.path.join(job_result_dir,f)
            #  logging.debug(p)
            # logging.debug(f'Process terminated. {job_title + " " + f.split(".")[-2]} file is missing.')
            # raise FileNotFoundError
    df_merged = reduce(lambda left, right: pd.merge(left, right, on=['pdb_position'],
                                                    how='outer'), dfs)
    contents = []
    for f in os.listdir(job_result_dir):
        try:
            if f.split(".")[-2] in ['disorder', 'binding']:
                p = os.path.join(job_result_dir, f)
                o = open(p, 'r')
                content = o.read().split('\n')
                o.close()
                contents.append(content)
        except:
            continue
    logging.debug(df_merged)

    df_merged.loc[df_merged['pdb_position'].isin(contents[0]), 'disorder'] = 'Disordered'
    df_merged.loc[df_merged['pdb_position'].isin(contents[1]), 'binding'] = 'Binding'
    df_merged['disorder'] = df_merged['disorder'].fillna('Non disordered')
    df_merged['binding'] = df_merged['binding'].fillna('Non binding')

    catalytic_sites_table = pd.read_csv(catalytic_sites)

    logging.debug(f'glycosylated: {0}')
    catalytic_sites_table['pdb_position'] = catalytic_sites_table['CHAIN ID'].astype(str) + catalytic_sites_table[
        'RESIDUE NUMBER'].astype(str) + "_" \
                                            + catalytic_sites_table['PDB ID'].str.upper()

    df_merged.loc[df_merged['pdb_position'].isin(catalytic_sites_table['pdb_position'].tolist()), 'catalysis'] = 'Catalytic'
    df_merged['catalysis'] = df_merged['catalysis'].fillna('None catalytic')
    logging.debug(f"columns after first merging={df_merged.columns.tolist()}")

    # calc AA FEATURES
    L=[c for c in df_merged.columns.tolist() if "wcn" in c]


    pdb_aa_col_name=[c for c in df_merged.columns.tolist() if "pdb_aa" in c]
    if len(pdb_aa_col_name)>0:
        pdb_aa_col_name= pdb_aa_col_name[0]
    else:
        pdb_aa_col_name = 'pdb_aa'
        df_merged[pdb_aa_col_name] = np.nan
    #if len(L) < 1:
     #   pdb_aa_col_name = 'pdb_aa'
    #else:
    #    pdb_aa_col_name = 'pdb_aa_x'

    if 'structure' not in df_merged.columns.tolist():
         df_merged['structure'] = np.nan



    df_merged.loc[df_merged[pdb_aa_col_name].isin(list('IVL')), 'aa_group_5'] = 'Aliphatic'
    df_merged.loc[df_merged[pdb_aa_col_name].isin(list('FYWH')), 'aa_group_5'] = 'Aromatic'
    df_merged.loc[df_merged[pdb_aa_col_name].isin(list('KRDE')), 'aa_group_5'] = 'Charged'
    df_merged.loc[df_merged[pdb_aa_col_name].isin(list('GACS')), 'aa_group_5'] = 'Tiny'
    df_merged.loc[df_merged[pdb_aa_col_name].isin(list('TMQNP')), 'aa_group_5'] = 'Diverse'
    df_merged.loc[df_merged[pdb_aa_col_name].isin(list('AGTSNQDEHRKP')), 'aa_group_HP'] = 'Polar'
    df_merged.loc[df_merged[pdb_aa_col_name].isin(list('CMFILVWY')), 'aa_group_HP'] = 'Hydrophobic'
    neigh_df = df_merged[['pdb_position', *[f for f in df_merged.columns.tolist() if 'neighbor_pos_' in f]]]

    neigh_d = dict(zip(neigh_df['pdb_position'],
                       neigh_df[[f for f in df_merged.columns.tolist() if 'neighbor_pos_' in f]].values.tolist()))

    logging.debug(f'glycosylated: {0}')
    # if not consurf_output:
    curr_residues = df_merged[pdb_aa_col_name].unique().tolist()
    curr_ss = df_merged['structure'].unique().tolist()
    curr_groups = df_merged['aa_group_5'].unique().tolist()
    curr_groups_HP = df_merged['aa_group_HP'].unique().tolist()

    aa_dicts = position2property_dict(df_merged,pdb_aa_col_name)
    ss_dicts = position2property_dict(df_merged,'structure')
    aa_groups_dicts = position2property_dict(df_merged,'aa_group_5')
    aa_groups_HP_dicts = position2property_dict(df_merged,'aa_group_HP')

    catalytic_sites_map = position2property_dict(df_merged,'catalysis')
    disordered_sites_map = position2property_dict(df_merged,'disorder')
    binding_sites_map = position2property_dict(df_merged,'binding')
    wcn_sc_map = position2property_dict(df_merged,'wcn_sc')
    wcn_ca_map = position2property_dict(df_merged,'wcn_ca')
    rsa_map = position2property_dict(df_merged,'rsa')
    logging.debug(f'glycosylated: {282}')

    print(len(catalytic_sites_map))
    print(len(neigh_d))
    print(len(df_merged))
    print([f for f in list(catalytic_sites_map.keys()) if f not in df_merged['pdb_position'].tolist()])
    print([f for f in df_merged['pdb_position'].tolist() if f not in list(catalytic_sites_map.keys())])
    print(df_merged['pdb_position'].tolist())
    df_merged = df_merged[~df_merged['pdb_position'].str.split("_").str[0].str[-1].isin(list(string.ascii_letters))]
    print(df_merged['pdb_position'].tolist())
    logging.debug(len(df_merged))
    if 'chain_x'in df_merged.columns.tolist() and 'chain_y' in df_merged.columns.tolist():
        df_merged = df_merged.dropna(subset=['chain_x','chain_y'])
    logging.debug(df_merged)
    logging.debug(len(df_merged))
    test_del=[[catalytic_sites_map[k] for k in neigh_d[row] if catalytic_sites_map.get(k, 0) == 'Catalytic'] for row in neigh_d if row in df_merged.pdb_position.tolist()]
    logging.debug(len(test_del))
    df_merged = df_merged.assign(
        total_catalytic_neigh=[[catalytic_sites_map[k] for k in neigh_d[row] if catalytic_sites_map.get(k, 0) == 'Catalytic']  for  row in neigh_d if row in df_merged.pdb_position.tolist()])
    df_merged = df_merged.assign(
        total_disordered_neigh=[[disordered_sites_map[k] for k in neigh_d[row] if disordered_sites_map.get(k, 0) == 'Disordered']  for  row in  neigh_d if row in df_merged.pdb_position.tolist()])
    df_merged = df_merged.assign(
        total_site_neigh=[[binding_sites_map[k] for k in neigh_d[row] if binding_sites_map.get(k, 0) == 'Binding']  for row in neigh_d if row in df_merged.pdb_position.tolist()])
    df_merged = df_merged.assign(
        median_wcn_ca_neigh=[[wcn_ca_map[k] for k in neigh_d[row] if wcn_ca_map.get(k, 0)] for row in neigh_d if row in df_merged.pdb_position.tolist()])
    df_merged = df_merged.assign(
        median_wcn_sc_neigh=[[wcn_sc_map[k] for k in neigh_d[row] if wcn_sc_map.get(k, 0)] for row in neigh_d if row in df_merged.pdb_position.tolist()])
    df_merged = df_merged.assign(
        median_rsa_neigh=[[rsa_map[k] for k in neigh_d[row] if rsa_map.get(k, 0)] for row in neigh_d if row in df_merged.pdb_position.tolist()])

    logging.debug(f'glycosylated: {318}')

    df_merged['total_catalytic_neigh'] = df_merged['total_catalytic_neigh'].apply(len)
    df_merged['total_disordered_neigh'] = df_merged['total_disordered_neigh'].apply(len)
    df_merged['total_site_neigh'] = df_merged['total_site_neigh'].apply(len)
    df_merged['median_wcn_ca_neigh'] = df_merged['median_wcn_ca_neigh'].apply(np.median)
    df_merged['median_wcn_sc_neigh'] = df_merged['median_wcn_sc_neigh'].apply(np.median)
    df_merged['median_rsa_neigh'] = df_merged['median_rsa_neigh'].apply(np.median)

    logging.debug(f'df: {df_merged}')

    get_neigh_features(aa_dicts, curr_residues, df_merged, neigh_d)

    get_neigh_features(ss_dicts, curr_ss, df_merged, neigh_d)


    get_neigh_features(aa_groups_dicts, curr_groups, df_merged, neigh_d)

    get_neigh_features(aa_groups_HP_dicts, curr_groups_HP, df_merged, neigh_d)
    if len(chain_list) > 1 and not contacting_residues_table.isnull().values.any():

        df_merged.loc[df_merged['pdb_position'].isin(contacting_residues_table[contacting_residues_table[2] == 'C']['pdb_position'].tolist()), 'protein_interaction'] = 'Contact'
        df_merged.loc[df_merged['pdb_position'].isin(contacting_residues_table[contacting_residues_table[2] == 'I'][
                                                         'pdb_position'].tolist()), 'protein_interaction'] = 'Interface'
        df_merged['protein_interaction'] = df_merged['protein_interaction'].fillna('Non interface')

        protein_interaction_dicts = position2property_dict(df_merged, 'protein_interaction')
        logging.debug(f'total_contact_neigh: {[[k for k in neigh_d[row] if protein_interaction_dicts.get(k, 0) == "Contact"] for row in neigh_d]}')
        logging.debug(f'total_contact_neigh: {len([[k for k in neigh_d[row] if protein_interaction_dicts.get(k, 0) == "Contact"] for row in neigh_d])}')
        logging.debug(f'total_contact_neigh: {[[protein_interaction_dicts[k] for k in neigh_d[row] if protein_interaction_dicts.get(k, 0) == "Contact"] for row in neigh_d]}')
        logging.debug(f'total_contact_neigh: {len([[protein_interaction_dicts[k] for k in neigh_d[row] if protein_interaction_dicts.get(k, 0) == "Contact"] for row in neigh_d])}')
        try:
            if not df_merged.shape[0]==len([[k for k in neigh_d[row] if protein_interaction_dicts.get(k, 0) == "Contact"] for row in neigh_d]):
                df_merged = df_merged.dropna(subset=['chain_x'],axis=0)
            df_merged['pdb_aa'] = df_merged['pdb_aa_x'].fillna(df_merged['pdb_aa_y'])
        except:
            pass
        logging.debug(f'df_merged: {df_merged}')
        debug_it = [[k for k in neigh_d[row] if protein_interaction_dicts.get(k, 0) == "Contact"] for row in neigh_d]
        debug_it = [item for sublist in debug_it for item in sublist]
        logging.debug(f'debug: {len(debug_it)}')
        logging.debug(f'debug: {[p for p in debug_it if p not in df_merged.pdb_position.tolist()]}')

        df_merged = df_merged.assign(total_contact_neigh=[[protein_interaction_dicts[k] for k in neigh_d[row] if protein_interaction_dicts.get(k, 0) == 'Contact'] for row in neigh_d if row in df_merged.pdb_position.tolist()])
        df_merged = df_merged.assign(total_interface_neigh=[[protein_interaction_dicts[k] for k in neigh_d[row] if protein_interaction_dicts.get(k, 0) == 'Interface'] for row in neigh_d if row in df_merged.pdb_position.tolist()])
        df_merged['total_contact_neigh'] = df_merged['total_contact_neigh'].apply(len)
        df_merged['total_interface_neigh'] = df_merged['total_interface_neigh'].apply(len)
    else:
        df_merged['protein_interaction'] = np.nan
        df_merged['total_contact_neigh'] = np.nan
        df_merged['total_interface_neigh'] = np.nan

    df_merged.loc[df_merged['pdb_position'].isin(glycosylated_positions), 'glycosylation'] = 'Glycosylated'
    df_merged['glycosylation'] = df_merged['glycosylation'].fillna('None glycosylated')

    glycosylation_dicts = position2property_dict(df_merged, 'glycosylation')

    df_merged = df_merged.assign(
        total_glycosylated_neigh=[[glycosylation_dicts[k] for k in neigh_d[row] if glycosylation_dicts.get(k, 0) == 'Glycosylated'] for
                          row in neigh_d if row in df_merged.pdb_position.tolist()])


    df_merged['total_glycosylated_neigh'] = df_merged['total_glycosylated_neigh'].apply(len)
    # Write protein strcuture feature set
    df_merged.to_csv(os.path.join(job_result_dir,job_title+"_feature_set.csv"),index=False)
    # Return protein strcuture feature set
    logging.debug(f'glycosylated: {347}')

    return df_merged




