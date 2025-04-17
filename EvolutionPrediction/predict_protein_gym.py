
import click

import os

laptop = False if os.path.exists("/groups/pupko/natannag/natan_git") else True
if laptop:
    path2github = "../"
    path2scannet = "C:/Users/natan/Documents/ScanNet_dev/"
    # path2scannet = "C:/Users/natan/Documents/ScanNet/"
    # path2scannet = "C:/Users/natan/Documents/missense_pathogenecity/ScanNet-main/ScanNet/"""
    path2evolutionprediction = "C:/Users/natan/Documents/EvolutionPrediction/"
else:
    path2scannet = "/groups/pupko/natannag/natan_git/ScanNet_dev"
import sys
sys.path.append(path2scannet)
import os
from ScanNet_dev.preprocessing import PDBio
import pandas as pd
from ScanNet_dev import predict_features
from Bio.PDB import PDBIO, PDBParser
import predict_aa_variance
import numpy as np
import sys
from Bio.SubsMat import MatrixInfo as matlist

data = r"D:\university\projects\nathan\evorator\evorator4_aa_variety_pred.csv"


@click.command()
@click.argument('evorator_df_path', type=str)
@click.argument('scannet_df_path', type=str)
@click.argument('output_dir', type=str)
@click.option('--predict-dms','-pd', is_flag=True)
@click.option('--write-cmds','-wc', is_flag=True)
def main(evorator_df_path,scannet_df_path, output_dir,predict_dms,write_cmds):

    if laptop:
        meta_data = pd.read_csv('ProteinGym_reference_file_substitutions.csv')
        meta_data2 = pd.read_csv('proteingym_to_uniprot_ids.tsv', sep='\t')
    else:
        meta_data = pd.read_csv('/bioseq/evorator/EvolutionPrediction/ProteinGym_reference_file_substitutions.csv')
        meta_data2 = pd.read_csv('/bioseq/evorator/EvolutionPrediction/proteingym_to_uniprot_ids.tsv', sep='\t')
    huang = True
    os.makedirs(output_dir, exist_ok=True)
    i = 0
    results_report = 'DMS_log.log'
    log_file = open(results_report, "a")
    log_file.write("\n")
    dms_dfs = []
    with open('proteingym_evorator_cmds.cmds', 'w') as evorator_cmds:
        if write_cmds:
            for f in os.listdir('substitutions'):
                # if i==15:
                #     break
                # try:
                uniprot_id = meta_data[meta_data.DMS_filename == f]['UniProt_ID'].values[0]
                print(uniprot_id)
                log_file.write(f"{uniprot_id}")
                log_file.write("\n")
                uniprot_id2 = meta_data2[meta_data2.From == uniprot_id]['Entry'].values[0]
                print(uniprot_id2)
                log_file.write(f"{uniprot_id2}")
                log_file.write("\n")

                path, chain = PDBio.getPDB(uniprot_id2)
                print(path)
                log_file.write(path)
                log_file.write("\n")
                if path is None:
                    continue
                OUTPUT_NAME = os.path.split(path)[-1].split('.')[0] + "_scannet_features.csv"
                OUTPATH = os.path.join(output_dir, OUTPUT_NAME)

                if os.path.exists(OUTPATH):
                    parser = PDBParser(QUIET=True)
                    if write_cmds:
                        with open(path, 'r') as pdb_file:
                            structure = parser.get_structure('struct', pdb_file)

                            chain_list = [c for c in structure.get_chains()]
                        print(path, chain, chain_list)
                        log_file.write(f"{path}, {chain}, {chain_list}")
                        log_file.write("\n")
                        dms_data = pd.read_csv(os.path.join('substitutions', f))
                        dms_seq = dms_data.mutant.str[0] + dms_data.mutated_sequence.str[1:]
                        print(uniprot_id2, dms_seq.values[0])
                        log_file.write(f"{uniprot_id2}, {dms_seq.values[0]}")
                        log_file.write("\n")

                        cmd = f'module load python/python-anaconda3.7-itaym;source activate /groups/pupko/natannag/conda/envs/NatanEnv;python /groups/pupko/natannag/consurf_n2v/huang/evorator_final_backup_270722.py "" "/groups/pupko/natannag/natan_git/EvolutionPrediction/proteingym_evorator/PDB/AF-{uniprot_id2}.pdb" {chain_list[0].id} /groups/pupko/natannag/consurf_n2v/huang/catalytic_sites.csv /groups/pupko/natannag/natan_git/EvolutionPrediction/proteingym_evorator/evorator_output/AF-{uniprot_id2} --orphan-prediction="True" --trained-regressor=/groups/pupko/natannag/consurf_n2v/huang/results_4_webserver/model_selection/svr_fs/EvoRator_SVR_final_model_PERFO_281121.joblib --consurfdb-query="" --consurf-output="" --prediction-task="" --job-title=AF-{uniprot_id2}\tevo_{uniprot_id2}\n'
                        print(cmd)
                        log_file.write(cmd)
                        log_file.write("\n")
                        evorator_cmds.write(cmd)
            else:
                print(uniprot_id, uniprot_id2)
                log_file.write(f"{uniprot_id}")
                log_file.write(f"{uniprot_id2}")
                log_file.write("\n")
                layers = [
                    'all_embedded_attributes_aa',  # dim 96 (small neighborhood)
                    'SCAN_filter_activity_aa',  # dim 128 (larger neighborhood ~13A)
                    'SCAN_filters_aa_embedded_1',  # dim 32 (low-dimensional projection from SCAN_filter_activity_aa)
                ]
                list_dictionary_features = predict_features.predict_features(path,
                                                                             model='ScanNet_PPI_noMSA',
                                                                             # PPBS model without evolution.
                                                                             layer=layers,
                                                                             # layers[0], # AA-scale spatio-chemical filters
                                                                             output_format='dictionary',
                                                                             biounit=True,
                                                                             permissive=False)
                print(list_dictionary_features.keys())
                log_file.write("\n")
                mapping_evorator_key = [k[1] + str(k[-1]) + "_" + os.path.split(path)[-1].split('.')[0] for k in
                                        list_dictionary_features.keys()]
                print(mapping_evorator_key[:10])
                log_file.write(f"{mapping_evorator_key[:10]}")
                log_file.write("\n")
                features_per_res = [np.concatenate(list_dictionary_features[k]) for k in list_dictionary_features]
                print(features_per_res[:10])
                log_file.write(f"{features_per_res[:10]}")
                log_file.write("\n")
                FINAL_DF = pd.DataFrame(features_per_res, index=mapping_evorator_key,
                                        columns=['scannet_feature' + str(i) for i in
                                                 range(len(features_per_res[0]))]).reset_index()
                FINAL_DF = FINAL_DF.rename(columns={'index': 'pdb_position'})
                FINAL_DF.to_csv(OUTPATH, index=False)
            i += 1

        elif predict_dms:
            evorator_df = pd.read_csv(evorator_df_path)
            evorator_df['pdb_position'] = evorator_df['pdb_position'].str.replace("-", "_")
            scannet_df = pd.read_csv(scannet_df_path)
            evo_sn_df = evorator_df.merge(scannet_df, on='pdb_position', indicator=True)
            print('shared predictions', evo_sn_df._merge.value_counts())
            log_file.write('shared predictions' + f'{evo_sn_df._merge.value_counts()}')
            log_file.write("\n")
            uniprot_id = scannet_df_path.split("/")[-1].split("_")[1]
            if laptop:
                evo_sn_df_path = os.path.join('proteingym_evorator', f'AF_{uniprot_id}_evorator_scannet_features.csv')
            else:
                evo_sn_df_path = os.path.join('/bioseq/evorator/EvolutionPrediction/proteingym_evorator/evorator_output', f'AF_{uniprot_id}_evorator_scannet_features.csv')
            if not os.path.exists(evo_sn_df_path):
                evo_sn_df.to_csv(evo_sn_df_path, index=False)
            print('PREDICTING')
            log_file.write('PREDICTING')
            log_file.write("\n")
            if laptop:
                if huang:
                    PRED_OUTPATH = predict_aa_variance.map_residue_variety_to_features([evo_sn_df_path,
                                                                                        "results/huang/ANN/trained_models/debug_cwv_evorator_only_ANN_n_hidden_147_trained_model_del.pkl",
                                                                                        "results/huang/ANN/trained_models/debug_cwv_evorator_only_ANN_n_hidden_147_fitted_preprocessor_del.pkl",
                                                                                        "200",
                                                                                        f'-o {uniprot_id}_evorator_scannet_huang'],
                                                                                       standalone_mode=False)
                else:
                    PRED_OUTPATH = predict_aa_variance.map_residue_variety_to_features([evo_sn_df_path,
                                                                                        "results/model_selection/deep_learning_for_proteinnet/trained_models/debug_to_4e5_data_pp_ppi_trained_model_del.h5",
                                                                                        "results/huang/ANN/trained_models/debug_cwv_ANN_n_hidden_147_fitted_preprocessor_del.pkl",
                                                                                        f'-o',f'{uniprot_id}_evorator_scannet_huang'],
                                                                                       standalone_mode=False)

            else:
                # PRED_OUTPATH = predict_aa_variance.map_residue_variety_to_features([evo_sn_df_path,
                #                                                                     "/bioseq/evorator/ml_results/model_selection/deep_learning_for_proteinnet/trained_models/debug_to_1e5_data_pp_ppi_trained_model_del.h5",
                #                                                                     "/bioseq/evorator/ml_results/model_selection/deep_learning_for_proteinnet/trained_models/debug_to_1e5_data_pp_ppi_fitted_preprocessor_del.pkl",
                #                                                                        '467',
                #                                                                     f'-o',f'{uniprot_id}_evorator_scannet_huang'],
                #                                                                    standalone_mode=False)

                PRED_OUTPATH = predict_aa_variance.map_residue_variety_to_features([evo_sn_df_path,
                                                                                    "/bioseq/evorator/ml_results/model_selection/deep_learning_for_proteinnet/trained_models/debug_to_all_data_pp_ppi_trained_model_del_2",
                                                                                    "/bioseq/evorator/ml_results/model_selection/deep_learning_for_proteinnet/trained_models/debug_to_all_data_pp_ppi_fitted_preprocessor_del.pkl",
                                                                                       '515',
                                                                                    f'-o',
                                                                                    f'{uniprot_id}_evorator_scannet_huang'],
                                                                                   standalone_mode=False)

            print(PRED_OUTPATH)
            log_file.write(PRED_OUTPATH)
            log_file.write("\n")
            evo_sn_df = pd.read_csv(PRED_OUTPATH)
            # print(evo_sn_df.columns.tolist())
            # print(evo_sn_df.columns.tolist())
            log_file.write("\n")

            # print(uniprot_id)
            # print(uniprot_id)
            log_file.write("\n")
            dms_id = meta_data2[meta_data2.Entry == uniprot_id]['From'].values[0]
            # print(dms_id)
            # print(dms_id)
            log_file.write("\n")
            # print(meta_data.UniProt_ID.value_counts())
            # print(meta_data.UniProt_ID.value_counts())
            log_file.write("\n")
            counts = meta_data.UniProt_ID.value_counts()
            # if counts[dms_id] == 1:
            print('easy - single data set')
            log_file.write('easy - single data set')
            log_file.write("\n")
            print(dms_id)
            log_file.write(f"{dms_id}")
            log_file.write("\n")
            dms_df_names = meta_data[meta_data.UniProt_ID == dms_id]['DMS_filename'].values
            print(dms_df_names)
            log_file.write(f"{dms_df_names}")
            log_file.write("\n")
            for dms_df_name in dms_df_names:
                # if os.path.exists('substitutions/' + dms_df_name.replace(".csv", "_evorator_scannet.csv")):
                #
                #     dms_df = pd.read_csv('substitutions/' + dms_df_name.replace(".csv", "_evorator_scannet.csv"))
                #     return dms_df
                if not laptop:
                    dms_df = pd.read_csv('/bioseq/evorator/EvolutionPrediction/substitutions/' + dms_df_name)
                else:
                    dms_df = pd.read_csv('substitutions/' + dms_df_name)
                print('Matching DMS DATA and EVORATOR_SCANNET output')
                log_file.write('Matching DMS DATA and EVORATOR_SCANNET output')
                log_file.write("\n")

                assert meta_data[meta_data.UniProt_ID == dms_id]['seq_len'].values[0] == evo_sn_df.shape[
                    0], f'dms sequence length {meta_data[meta_data.UniProt_ID == dms_id]["seq_len"].values[0]} does not match evorator_scannet {evo_sn_df.shape[0]}'
                print('Identical sequence lengths')
                log_file.write('Identical sequence lengths')
                log_file.write("\n")
                # print(dms_df)
                # print(dms_df)
                log_file.write("\n")
                # if not meta_data[meta_data.UniProt_ID == dms_id]['includes_multiple_mutants'].values[0]:
                print('easy - single mutatiotns')
                log_file.write('easy - single mutatiotns')
                log_file.write("\n")
                evo_sn_df['map_key'] = evo_sn_df['index'].str.split("_").str[0]
                evo_sn_df['map_key'] = evo_sn_df['map_key'].str.extract('(\d+)').astype(int)
                evo_sn_df['map_key'] = evo_sn_df.pdb_aa + evo_sn_df['map_key'].astype(str)

                # print(evo_sn_df.map_key)
                # print(evo_sn_df.map_key)
                log_file.write("\n")
                # print(evo_sn_df.columns.tolist())
                # print(evo_sn_df.columns.tolist())
                log_file.write("\n")
                print('Matching DMS DATA and EVORATOR_SCANNET output')
                log_file.write('Matching DMS DATA and EVORATOR_SCANNET output')
                log_file.write("\n")
                mismatch = False
                if not meta_data[meta_data.UniProt_ID == dms_id]['target_seq'].values[0] == ''.join(
                        evo_sn_df.pdb_aa.tolist()):
                    mismatch = True
                    print("MISMATCH")
                    log_file.write("MISMATCH")
                    log_file.write("\n")
                    print(
                        f'dms sequence {meta_data[meta_data.UniProt_ID == dms_id]["target_seq"].values[0]} does not match evorator_scannet {"".join(evo_sn_df.pdb_aa.tolist())}')
                    log_file.write(
                        f'dms sequence {meta_data[meta_data.UniProt_ID == dms_id]["target_seq"].values[0]} does not match evorator_scannet {"".join(evo_sn_df.pdb_aa.tolist())}')
                    log_file.write("\n")
                if not mismatch:
                    print('Identical target sequence')
                    log_file.write('Identical target sequence')
                    log_file.write("\n")
                pred_dict = {}
                aa_list = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V',
                           'W', 'Y']
                # print(evo_sn_df[aa_list])
                # print(evo_sn_df[aa_list])
                log_file.write("\n")
                for k in evo_sn_df['map_key'].tolist():
                    for c in aa_list:
                        pred_dict[k + c] = float(evo_sn_df.loc[evo_sn_df.map_key == k, c])
                # S136C:N137D:V143I:E144N:I160V:P161T:F163V
                multi_mut_d = {}
                for mut_lst in dms_df['mutant'].str.split(":").tolist():
                    for mut in mut_lst:
                        if ":".join(mut_lst) not in multi_mut_d:
                            multi_mut_d[":".join(mut_lst)] = [pred_dict.get(mut, 1)]
                        else:
                            multi_mut_d[":".join(mut_lst)].append(pred_dict.get(mut, 1))
                multi_mut_multiplied_d = {}
                for k in multi_mut_d:
                    tmp = 0
                    for score in multi_mut_d[k]:
                        tmp += score
                    multi_mut_multiplied_d[k] = tmp
                blosum_matrix = matlist.blosum62

                print(list(multi_mut_multiplied_d.items())[:10])
                log_file.write(f"{list(multi_mut_multiplied_d.items())[:10]}")
                log_file.write("\n")
                blosum_pred = {dms_df.loc[i, 'mutant']: blosum_matrix.get(
                    (dms_df.loc[i, 'mutant'][0], dms_df.loc[i, 'mutant'][-1]),
                    blosum_matrix.get((dms_df.loc[i, 'mutant'][-1], dms_df.loc[i, 'mutant'][0]))) for i in
                    range(dms_df.shape[0])}
                print(list(blosum_pred.items())[:10])
                log_file.write(f"{list(blosum_pred.items())[:10]}")
                log_file.write("\n")
                dms_df['evorator_pred'] = dms_df['mutant'].map(multi_mut_multiplied_d)
                dms_df['blosum_pred'] = dms_df['mutant'].map(blosum_pred)

                if mismatch:
                    dms_df['evorator_pred'] = dms_df['evorator_pred'].fillna(191189)
                    dms_df['blosum_pred'] = dms_df['blosum_pred'].fillna(191189)
                    n_match = dms_df[dms_df.evorator_pred != 191189].shape[0]
                    n_mismatch = dms_df[dms_df.evorator_pred == 191189].shape[0]
                    n_total = dms_df.shape[0]
                    print(" total positions=", n_total)
                    log_file.write(f" total positions=, {n_total}")
                    print(" mismatch position=", n_mismatch)
                    log_file.write(f" mismatch position=, {n_mismatch}")
                    print(f"% mismatch position= {round(100 * (n_mismatch / n_total), 2)}")
                    log_file.write(f"% mismatch position= {round(100 * (n_mismatch / n_total), 2)}")
                    log_file.write("\n")
                    dms_df = dms_df[dms_df.evorator_pred != 191189]
                    log_file.write("\n")
                    n_total = dms_df.shape[0]

                # print(pred_dict.keys(), dms_df['mutant'])
                # print(pred_dict.keys(), dms_df['mutant'])
                log_file.write("\n")
                # print(dms_df[['mutant', 'evorator_pred']].head(3000))
                # print(dms_df[['mutant', 'evorator_pred']].head(3000))
                log_file.write("\n")
                # dms_df[RESPONSE] = dms_df[RESPONSE].apply(np.log)
                # dms_df['evorator_pred'] = dms_df['evorator_pred'].apply(np.log)
                print(dms_df[['evorator_pred', 'blosum_pred',
                              # log_file.write(dms_df[['evorator_pred', 'blosum_pred',
                              #               log_file.write("\n")
                              *'DMS_score_bin,Tranception_L_no_retrieval,Tranception_S_retrieval,Tranception_M_retrieval,Tranception_L_retrieval,EVE_single,EVE_ensemble,MSA_Transformer_single,MSA_Transformer_ensemble,ESM1v_single,ESM1v_ensemble,Wavenet,DeepSequence_single,DeepSequence_ensemble,Site_Independent,EVmutation,RITA_s,RITA_m,RITA_l,RITA_xl,RITA_ensemble,Progen2_small,Progen2_medium,Progen2_base,Progen2_large,Progen2_xlarge,Progen2_ensemble,Ensemble_Tranception_EVE'.split(
                                  ','), 'DMS_score']].dropna().corr(method='spearman').iloc[:, -1])
                log_file.write(
                    f"{dms_df[['evorator_pred', 'blosum_pred', *'DMS_score_bin,Tranception_L_no_retrieval,Tranception_S_retrieval,Tranception_M_retrieval,Tranception_L_retrieval,EVE_single,EVE_ensemble,MSA_Transformer_single,MSA_Transformer_ensemble,ESM1v_single,ESM1v_ensemble,Wavenet,DeepSequence_single,DeepSequence_ensemble,Site_Independent,EVmutation,RITA_s,RITA_m,RITA_l,RITA_xl,RITA_ensemble,Progen2_small,Progen2_medium,Progen2_base,Progen2_large,Progen2_xlarge,Progen2_ensemble,Ensemble_Tranception_EVE'.split(','), 'DMS_score']].dropna().corr(method='spearman').iloc[:, -1]}")
                X = dms_df[['evorator_pred', 'blosum_pred',
                            *'Tranception_L_no_retrieval,Tranception_S_retrieval,Tranception_M_retrieval,'
                             'Tranception_L_retrieval,EVE_single,EVE_ensemble,MSA_Transformer_single,MSA_Transformer_ensemble,ESM1v_single'
                             ',ESM1v_ensemble,Wavenet,DeepSequence_single,DeepSequence_ensemble,Site_Independent,EVmutation,RITA_s,RITA_m,'
                             'RITA_l,RITA_xl,RITA_ensemble,Progen2_small,Progen2_medium,Progen2_base,'
                             'Progen2_large,Progen2_xlarge,Progen2_ensemble,Ensemble_Tranception_EVE'.split(',')]]

                y = dms_df.DMS_score
                # import statsmodels.formula.api as sm
                # res1 = sm.ols(formula="DMS_score ~ evorator_pred", data=dms_df).fit()
                # print('evorator_scannet',res1.rsquared_adj)
                # print('evorator_scannet',res1.rsquared_adj)
                log_file.write("\n")
                # res2 = sm.ols(formula="DMS_score ~ Ensemble_Tranception_EVE", data=dms_df).fit()
                # print('tranception_eve',res2.rsquared_adj)
                # print('tranception_eve',res2.rsquared_adj)
                log_file.write("\n")
                # res12 = sm.ols(formula="DMS_score ~ evorator_pred+Ensemble_Tranception_EVE", data=dms_df).fit()
                # print('evorator_scannet_tranception_eve',res12.rsquared_adj)
                # print('evorator_scannet_tranception_eve',res12.rsquared_adj)
                log_file.write("\n")
                #
                dms_df['uniprot_id'] = uniprot_id
                dms_df['dms_id'] = dms_id
                # dms_df['evorator_scannet_r2_adj'] =  res1.rsquared_adj
                # dms_df['tranception_eve_r2_adj'] =  res2.rsquared_adj
                # dms_df['evorator_scannet_tranception_eve_r2_adj'] =  res12.rsquared_adj
                # MSA_num_seqs = meta_data[meta_data.UniProt_ID == dms_id]['MSA_num_seqs'].values[0]
                # with open(results_report,'a') as ff:
                #     ff.write(uniprot_id+'\n')
                #     ff.write('correlation'+'\n')
                #     ff.write('evorator_pred  Tranception_L_no_retrieval Tranception_S_retrieval Tranception_M_retrieval Tranception_L_retrieval EVE_single EVE_ensemble MSA_Transformer_single MSA_Transformer_ensemble ESM1v_single ESM1v_ensemble Wavenet DeepSequence_single DeepSequence_ensemble Site_Independent EVmutation RITA_s RITA_m RITA_l RITA_xl RITA_ensemble Progen2_small Progen2_medium Progen2_base Progen2_large Progen2_xlarge Progen2_ensemble Ensemble_Tranception_EVE'+'\n')
                #     ff.write(' '.join(dms_df[['evorator_pred',
                #               *'Tranception_L_no_retrieval,Tranception_S_retrieval,Tranception_M_retrieval,Tranception_L_retrieval,EVE_single,EVE_ensemble,MSA_Transformer_single,MSA_Transformer_ensemble,ESM1v_single,ESM1v_ensemble,Wavenet,DeepSequence_single,DeepSequence_ensemble,Site_Independent,EVmutation,RITA_s,RITA_m,RITA_l,RITA_xl,RITA_ensemble,Progen2_small,Progen2_medium,Progen2_base,Progen2_large,Progen2_xlarge,Progen2_ensemble,Ensemble_Tranception_EVE'.split(
                #                   ','), 'DMS_score']].dropna().corr(method='spearman').iloc[:,-1].astype(str).values.tolist())+'\n')
                #     ff.write('uniprot_id | dms_id | evorator_scannet | tranception_eve | evorator_scannet_tranception_eve | MSA_perc_cov'+'\n')
                #     ff.write(uniprot_id + " | " + dms_id + " | "+ str(res1.rsquared_adj)+' | '+str(res2.rsquared_adj)+' | '+str(res12.rsquared_adj)+ " | "+str(MSA_num_seqs)+'\n')
                print('writing_results')
                log_file.write('writing_results')
                log_file.write("\n")
                if not laptop:
                    dms_df.to_csv('/bioseq/evorator/EvolutionPrediction/substitutions/' + dms_df_name.replace(".csv","_evorator_scannet.csv"),index=False)
                else:
                    dms_df.to_csv('substitutions/' + dms_df_name.replace(".csv","_evorator_scannet.csv"),index=False)
                print('finsihed writing')
                log_file.write('finsihed writing')
                log_file.write("\n")
                dms_dfs.append(dms_df)
            log_file.close()
            return dms_dfs
                # else:
            #     print('hard - multiple mutatiotns')
            #     print('hard - multiple mutatiotns')
            #     pass

            # except Exception as e:
            #     print(e)
            #     print(e)

        # return uniprot id to best aligned pdb id
        #
        #
        # predict proteingym records that do not share CATH with HUANG


if __name__ == '__main__':
    import os
    laptop = False if os.path.exists("/groups/pupko/natannag/natan_git") else True
    if laptop:
        path2github = "../"
        path2scannet = "C:/Users/natan/Documents/ScanNet_dev/"
        # path2scannet = "C:/Users/natan/Documents/ScanNet/"
        # path2scannet = "C:/Users/natan/Documents/missense_pathogenecity/ScanNet-main/ScanNet/"""
        path2evolutionprediction = "C:/Users/natan/Documents/EvolutionPrediction/"
    sys.path.append(path2scannet)
    import os
    from ScanNet_dev.preprocessing import PDBio
    import pandas as pd
    from ScanNet_dev import predict_features
    from Bio.PDB import PDBIO,PDBParser

    main()
# DMS natan_git/dataset to uniprot id
