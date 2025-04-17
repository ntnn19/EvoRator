import predict_protein_gym
import os
import pandas as pd
import logging
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


logging.basicConfig(
        filename=f'evorator_scannet_DMS_211122.log',
        level=logging.DEBUG,
        format='%(asctime)s %(levelname)-8s %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S')
if laptop:
    meta_data2 = pd.read_csv('proteingym_to_uniprot_ids.tsv', sep='\t')
    meta_data = pd.read_csv('ProteinGym_reference_file_substitutions.csv')
else:
    meta_data2 = pd.read_csv('/bioseq/evorator/EvolutionPrediction/proteingym_to_uniprot_ids.tsv', sep='\t')
    meta_data = pd.read_csv('/bioseq/evorator/EvolutionPrediction/ProteinGym_reference_file_substitutions.csv')

uniprot_ids = meta_data2['Entry'].unique()
dms_dfs_train = []
dms_dfs_test = []
i = 0
# for q in uniprot_ids:
dms_uniprot_ids_multi_datasets = (meta_data.UniProt_ID.value_counts()[meta_data.UniProt_ID.value_counts() > 1]).index
# print(meta_data.UniProt_ID.value_counts()[meta_data.UniProt_ID.value_counts()>1])
# print(meta_data[meta_data.UniProt_ID.value_counts()>1])
dms_uniprot_ids_multi = meta_data[meta_data.includes_multiple_mutants == True].UniProt_ID
uniprot_ids_multi = meta_data2[meta_data2['From'].isin(dms_uniprot_ids_multi)]['Entry'].unique()
uniprot_ids_multi_datasets = meta_data2[meta_data2['From'].isin(dms_uniprot_ids_multi_datasets)]['Entry'].unique()
uniprot_ids_single = meta_data2[~meta_data2['From'].isin(dms_uniprot_ids_multi)]['Entry'].unique()
print('single', uniprot_ids_single.shape)
print('multi', uniprot_ids_multi.shape)
print('multi data', uniprot_ids_multi_datasets.shape)
print('multi data', uniprot_ids_multi_datasets)
print('multi data', meta_data2)
for q in dms_uniprot_ids_multi_datasets:  # debug multiple mutations
    print(q)
    # if i==1:
    #     break
    try:
        if i % 5 == 0:
            if not laptop:
                dms_dfs_test.append(predict_protein_gym.main([
                                                                 f'/bioseq/evorator/EvolutionPrediction/proteingym_evorator/evorator_output/AF-{q}_A/AF_{q}A_features_and_predictions.csv',
                                                                 f'/bioseq/evorator/EvolutionPrediction/proteingym_scannet/AF_{q}_scannet_features.csv',
                                                                 'proteingym_evorator_scannet', '-pd'],
                                                             standalone_mode=False))

            else:
                dms_dfs_test.append(predict_protein_gym.main(
                    [f'proteingym_evorator/evorator_output/AF-{q}/AF-{q}_features_and_predictions.csv',
                     f'proteingym_scannet/AF_{q}_scannet_features.csv',
                     'proteingym_evorator_scannet', '-pd'], standalone_mode=False))

        else:

            if not laptop:
                dms_dfs_train.append(predict_protein_gym.main(
                [f'/bioseq/evorator/EvolutionPrediction/proteingym_evorator/evorator_output/AF-{q}_A/AF_{q}A_features_and_predictions.csv',
                 f'/bioseq/evorator/EvolutionPrediction/proteingym_scannet/AF_{q}_scannet_features.csv',
                 'proteingym_evorator_scannet', '-pd'], standalone_mode=False))

            else:
                dms_dfs_train.append(predict_protein_gym.main(
                    [f'proteingym_evorator/evorator_output/AF-{q}/AF-{q}_features_and_predictions.csv',
                     f'proteingym_scannet/AF_{q}_scannet_features.csv',
                     'proteingym_evorator_scannet', '-pd'], standalone_mode=False))
        # i+=1
    except Exception as e:
        print('ERROR')
        logging.debug('ERROR')
        print(e)
        logging.debug(e)

for q in [*uniprot_ids_multi, *uniprot_ids_single]:  # debug multiple mutations
    print(q)
    # if i==20:
    #     break
    try:
        if i % 5 == 0:
            if not laptop:

                dms_dfs_test.append(predict_protein_gym.main([
                                                                 f'/bioseq/evorator/EvolutionPrediction/proteingym_evorator/evorator_output/AF_{q}_A/AF-{q}A_features_and_predictions.csv',
                                                                 f'/bioseq/evorator/EvolutionPrediction/proteingym_scannet/AF_{q}_scannet_features.csv',
                                                                 'proteingym_evorator_scannet', '-pd'],
                                                             standalone_mode=False))

            else:
                dms_dfs_test.append(predict_protein_gym.main(
                    # [f'proteingym_evorator/old_evorator_output/AF-{q}/AF-{q}_features_and_predictions.csv',

                    [f'proteingym_evorator/old_evorator_output/AF-{q}_features_and_predictions.csv',
                     f'proteingym_scannet/AF_{q}_scannet_features.csv',
                     'proteingym_evorator_scannet', '-pd'], standalone_mode=False))

        else:

            if not laptop:
                dms_dfs_train.append(predict_protein_gym.main(
                    [
                        f'/bioseq/evorator/EvolutionPrediction/proteingym_evorator/evorator_output/AF_{q}_A/AF-{q}A_features_and_predictions.csv',
                        f'/bioseq/evorator/EvolutionPrediction/proteingym_scannet/AF_{q}_scannet_features.csv',
                        'proteingym_evorator_scannet', '-pd'], standalone_mode=False))

            else:
                dms_dfs_train.append(predict_protein_gym.main(
                    # [f'proteingym_evorator/old_evorator_output/AF-{q}/AF-{q}_features_and_predictions.csv',
                    [f'proteingym_evorator/old_evorator_output/AF-{q}_features_and_predictions.csv',
                     f'proteingym_scannet/AF_{q}_scannet_features.csv',
                     'proteingym_evorator_scannet', '-pd'], standalone_mode=False))
        i += 1
    except Exception as e:
        print('ERROR')
        logging.debug('ERROR')
        print(e)
        logging.debug(e)

# print(len(dms_dfs_train),len(dms_dfs_test))
# print(dms_dfs_train[:2],dms_dfs_test[:2])
dms_dfs_train_flat = []
for df in dms_dfs_train:
    if type(df) is list:
        for d in df:
            dms_dfs_train_flat.append(d)
    else:
        dms_dfs_train_flat.append(df)

dms_dfs_test_flat = []
for df in dms_dfs_test:
    if type(df) is list:
        for d in df:
            dms_dfs_test_flat.append(d)
    else:
        dms_dfs_test_flat.append(df)

dms_df_fina = pd.concat([*dms_dfs_train_flat,*dms_dfs_test_flat],axis=0)
if laptop:
    dms_df_fina.to_csv("FINAL_DMS_RESULTS_OLD_EVORATOR.csv")
else:
    dms_df_fina.to_csv("/bioseq/evorator/EvolutionPrediction/FINAL_DMS_RESULTS_NEW_EVORATOR.csv")
# print(dms_df_train[['Ensemble_Tranception_EVE','evorator_pred','DMS_score']].corr(method='spearman'))
# dms_df_test = pd.concat(dms_dfs_test,axis=0)
# X_train , y_train= dms_df_train[['evorator_pred']], dms_df_train['DMS_score']
# X_test , y_test= dms_df_test[['evorator_pred']], dms_df_test['DMS_score']
# from sklearn.linear_model import LinearRegression
# from sklearn.svm import SVR
# from sklearn.preprocessing import StandardScaler
# from sklearn.metrics import r2_score
# lr = LinearRegression()
# svr = SVR()
# ss = StandardScaler()
# ss.fit(X_train)
# X_train=  ss.transform(X_train)
# X_test= ss.transform(X_test)
# lr.fit(X_train,y_train)
# print('training score',lr.score(X_train,y_train))
# y_hat = lr.predict(X_test)
# print('test score',r2_score(y_test,y_hat))
#
# svr.fit(X_train,y_train)
# print('training score',svr.score(X_train,y_train))
# y_hat = svr.predict(X_test)
# print('test score',r2_score(y_test,y_hat))
