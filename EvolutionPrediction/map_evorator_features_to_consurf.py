import click
@click.command()
@click.argument('list_of_consurfdb_queries', type=click.Path(exists=False)) # pairs of <pdb_id> , <chain_od>
@click.argument('outdir', type=click.Path(exists=None))
@click.argument('n', type=int) # number of queries to sample from consurf-db
@click.argument('m', type=int) # number of queries to sample from consurf-db
def main(list_of_consurfdb_queries,outdir,n,m):
    assert (m > n and  m < 127024) and n>=0, 'number of queries must be 1 <= n < 127,024'
    os.makedirs(outdir,exist_ok=True)
    # obtain conservation profiles and pdb files from consurf-db
    input = open(list_of_consurfdb_queries, 'r')
    content = input.readlines()
    input.close()
    import random
    # random_queries = random.sample(range(len(content)), n) #list(range(len(content)) curr implementation avoids calculating the number of lines
    content = list(np.array(content)[n:m])
    print(content[:10])
    AA = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

    # content = list(np.array(content)[[67514]])  # debug cif
    # content = list(np.array(content)[[29620]])  # debug multichain
    # content = list(np.array(content)[[?????]])  # debug single chain
    for line in content:
        line = line.strip().split()
        unique_consurfdb_query, path = line[0], line[1]
        chain = unique_consurfdb_query[4:]
        job_title = unique_consurfdb_query[:4]
        print("Obtaining conservation profile and coordinate file from ConSurf-DB for:", unique_consurfdb_query)
        tmp_outdir = os.path.join(outdir, unique_consurfdb_query)
        consurf_db_file = os.path.join(path, unique_consurfdb_query + '_consurf_summary.txt')
        consurf_pssm_file = os.path.join(path, unique_consurfdb_query + '_msa_positional_aa_frequency.csv')
        FINAL_OUTPUT = os.path.join(tmp_outdir, unique_consurfdb_query + '_evorator_scannet_pssm_consurf.csv')
        pdb_file = os.path.join(path, job_title + '.pdb')
        cif_file = os.path.join(path, job_title + '.cif')
        os.makedirs(tmp_outdir, exist_ok=True)
        if not os.path.exists(FINAL_OUTPUT):
            try:
                if not os.path.exists(os.path.join(tmp_outdir, unique_consurfdb_query + '_consurf_summary.txt')):
                    subprocess.check_output(f'ssh bioseq@powerweb1 cp {consurf_db_file} {tmp_outdir}', shell=True)
                if not os.path.exists(os.path.join(tmp_outdir, unique_consurfdb_query + '_msa_positional_aa_frequency.csv')):
                    subprocess.check_output(f'ssh bioseq@powerweb1 cp {consurf_pssm_file} {tmp_outdir}', shell=True)
            except Exception as e:
                print(e)
                continue

            try:
                if not os.path.exists(os.path.join(tmp_outdir, job_title + '.pdb')):
                    subprocess.check_output(f'ssh bioseq@powerweb1 cp {pdb_file} {tmp_outdir}', shell=True)
            except Exception as e:
                print(e)
                try:
                    subprocess.check_output(f'ssh bioseq@powerweb1 cp {cif_file} {tmp_outdir}', shell=True)
                except Exception as e:
                    print(e)
                    continue


            # EvoRator feature extraction
            # maybe implement chain extraction before feature extraction
            cmd = f'module load python/python-anaconda3.7-itaym;source activate /groups/pupko/natannag/conda/envs/NatanEnv;python /bioseq/evorator/EvoRator/evorator_final_backup_270722.py "" {os.path.join(tmp_outdir,job_title + ".xxx")} {chain} /groups/pupko/natannag/natan_git/EvoRator/data/catalytic_sites.csv {tmp_outdir} --orphan-prediction="True" --trained-regressor=/groups/pupko/natannag/natan_git/EvoRator/evorator_model/EvoRator_SVR_final_model_PERFO_281121.joblib --consurfdb-query="" --consurf-output="" --prediction-task="" --job-title={job_title}'
            flg = True
            if flg:
                print("Extracting EvoRator features for:", unique_consurfdb_query)
                if os.path.exists(os.path.join(tmp_outdir,job_title + ".pdb")):

                    try:
                        # cmd = f'module load python/python-anaconda3.7-itaym;source activate /groups/pupko/natannag/conda/envs/NatanEnv;python /bioseq/evorator/EvoRator/evorator_final_backup_270722.py "" {os.path.join(tmp_outdir,job_title + ".pdb")} {chain} /groups/pupko/natannag/natan_git/EvoRator/data/catalytic_sites.csv {tmp_outdir} --orphan-prediction="True" --trained-regressor=/groups/pupko/natannag/natan_git/EvoRator/evorator_model/EvoRator_SVR_final_model_PERFO_281121.joblib --consurfdb-query="" --consurf-output="" --prediction-task="" --job-title={job_title}'
                        cmd = cmd.replace('xxx','pdb')
                        print(cmd)
                        subprocess.check_output(cmd, shell=True)
                    except Exception as e:
                        print(e)
                        continue

                elif os.path.exists(os.path.join(tmp_outdir,job_title + ".cif")):
                    try:
                        # cmd = f'module load python/python-anaconda3.7-itaym;source activate /groups/pupko/natannag/conda/envs/NatanEnv;python /bioseq/evorator/EvoRator/evorator_final_backup_270722.py "" {os.path.join(tmp_outdir,job_title + ".pdb")} {chain} /groups/pupko/natannag/natan_git/EvoRator/data/catalytic_sites.csv {tmp_outdir} --orphan-prediction="True" --trained-regressor=/groups/pupko/natannag/natan_git/EvoRator/evorator_model/EvoRator_SVR_final_model_PERFO_281121.joblib --consurfdb-query="" --consurf-output="" --prediction-task="" --job-title={job_title}'
                        cmd = cmd.replace('xxx', 'cif')
                        print(cmd)
                        subprocess.check_output(cmd, shell=True)

                    except Exception as e:
                        print(e)
                        continue

            # Scannet feature extraction
            flg = True
            if flg:
                print("Extracting ScanNet features for:", unique_consurfdb_query)
                try:
                    cmd=f'module unload python/python-anaconda3.7-itaym; module load miniconda/miniconda3-4.7.12-environmentally; conda activate /groups/pupko/natannag/conda-envs/py_scannet; python /bioseq/evorator/ScanNet_dev/feature_extraction_scannet.py {tmp_outdir} {tmp_outdir}'
                    print(cmd)
                    subprocess.check_output(cmd,shell=True)

                except Exception as e:
                    print(e)
                    continue

            # map ScanNet features and conservation profiles to residues analyzed by EvoRator
            try:
                with open(os.path.join(tmp_outdir,'evorator_output_path.txt'),'r') as f1:
                    evorator_df_path = f1.read().strip()
                    evorator_df = pd.read_csv(evorator_df_path)
                    evorator_df = evorator_df.rename(columns={'predicted_score': 'evorator_conservation_score',
                                                              'total_nan_neigh': 'total_non_interface_non_contact_neigh',
                                                              'total_ss_nan_neigh': 'total_ss_loop_neigh'})
                    for c in evorator_df.columns.tolist():
                        if 'chain' in c and c != 'chain':
                            evorator_df["chain"] = evorator_df['chain'].fillna(evorator_df[c])
                        if 'pdb_aa' in c and c != 'pdb_aa':
                            evorator_df["pdb_aa"] = evorator_df['pdb_aa'].fillna(evorator_df[c])
                    columns_not_used_in_analysis = ['Unnamed: 0', 'index', 'pdb_id', 'pdb', 'pdb_x', 'pdb_y', 'chain_x',
                                                    'chain_y',
                                                    'pdb_aa_y', 'pdb_aa_x', 'normalized_score', 'ResInd', 'Score',
                                                    'protein_interaction', 'total_non_interface_non_contact_neigh',
                                                    'total_contact_neigh', 'total_non_interface_neigh'
                                                                           'zr4s_JC'] \
                                                   + [f for f in evorator_df.columns.tolist()
                                                      if f.startswith('neighbor_pos_') or f.startswith("n2v")]


                with open(os.path.join(tmp_outdir,'scannet_output_path.txt'),'r') as f2:
                    scannet_df_path = f2.read().strip()
                    scannet_df = pd.read_csv(scannet_df_path)

                # Process conservation profiles
                conservation_df = pd.read_csv(consurf_db_file, sep='\t', skiprows=15, skipfooter=4,
                                              header=None).dropna(axis=1)
                pssm_df = pd.read_csv(consurf_pssm_file, skiprows=4)


                conservation_df.columns = ['POS', 'SEQ', '3LATOM', 'consurf_conservation_score', 'color',
                                           'confidence_interval',
                                           'confidence_colors',
                                           'msa_data', 'residue_variety']
                conservation_df['chain'] = conservation_df['3LATOM'].str.split(":").str[-1]
                conservation_df['pdb_position_consurf'] = conservation_df['3LATOM'].str.split(":").str[0].str.extract('(\d+)')
                conservation_df['pdb_position'] = conservation_df['chain'].astype(str) + conservation_df['pdb_position_consurf'].astype(str) + "_" + os.path.split(evorator_df_path)[-1].split("_")[0][:4]
                position_d = dict(zip(conservation_df['POS'], conservation_df['pdb_position']))
                pssm_df['pdb_position'] = pssm_df.pos.map(position_d)
                pssm_df[AA] = pssm_df[AA] / 100
                pssm_df[AA] = pssm_df[AA].div(pssm_df[AA].sum(axis=1), 'index')

                print(conservation_df.head())
                print(pssm_df.head())

                # Merge tables
                print('conservation df size =',conservation_df.shape)
                print('pssm df size =',pssm_df.shape)
                print('evorator df size =',evorator_df.shape)
                print('scannet df size =',scannet_df.shape)

                final_df = evorator_df.loc[:, ~evorator_df.columns.isin(columns_not_used_in_analysis)].merge(
                    pssm_df[['pdb_position', *AA]],
                    on='pdb_position', indicator=True, how='left')

                final_df = final_df.loc[:, ~final_df.columns.isin(['_merge'])].merge(
                    conservation_df[['pdb_position', 'consurf_conservation_score', 'confidence_interval']],
                    on='pdb_position', indicator=True, how='left')


                final_df = final_df.loc[:, ~final_df.columns.isin(['_merge'])].merge(scannet_df, on='pdb_position',indicator=True, how='left')


                final_df['split_by_CATH'] = os.path.split(evorator_df_path)[-1].split("_")[0][:4].lower() + "_" + final_df['chain']
                print(final_df.head())
                print('final df size =',final_df.shape)
                final_df.to_csv(FINAL_OUTPUT, index=False)

            except Exception as e:
                print(e)
                continue


if __name__ == '__main__':
    import sys
    import subprocess
    import os
    import pickle

    laptop = False if os.path.exists("/groups/pupko/natannag/natan_git") else True
    if laptop:
        path2github = "../"
        path2scannet = "C:/Users/natan/Documents/ScanNet_dev/"
        path2evolutionprediction = "C:/Users/natan/Documents/EvolutionPrediction/"
    else:
        path2github = "/groups/pupko/natannag/natan_git"
        path2scannet = '/groups/pupko/natannag/natan_git/ScanNet_dev/'
        path2evolutionprediction = '/groups/pupko/natannag/natan_git/EvolutionPrediction/'

    sys.path.append(path2github)
    sys.path.append(path2scannet)

    from ScanNet_dev import predict_features
    # sys.path.append('/bioseq/evorator/auxiliaries')
    import auxiliaries.evorator_CONSTANTS as CONSTS  # from /effectidor/auxiliaries
    import numpy as np
    import pandas as pd
    main()

