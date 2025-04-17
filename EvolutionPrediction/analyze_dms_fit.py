import pandas as pd
import datatable as dt
import reduce_memory_usage


# adding auc data to spearman
def get_final_df(df_merged_with_bins_for_analysis, c_flt_tmp):
    df_merged_with_bins_for_analysis = df_merged_with_bins_for_analysis.replace('auc1', 'evorator_pred')
    df_merged_with_bins_for_analysis = df_merged_with_bins_for_analysis.replace('auc2', 'Ensemble_Tranception_EVE')
    df_merged_with_bins_for_analysis = df_merged_with_bins_for_analysis.replace('auc12', 'm12_lrm_pred')
    auc_d = dict(zip(df_merged_with_bins_for_analysis['variable'], df_merged_with_bins_for_analysis['value']))
    auc_d = {}
    auc_dfs = []
    for c in df_merged_with_bins_for_analysis['DMS_id'].unique():
        tmp_df = df_merged_with_bins_for_analysis[df_merged_with_bins_for_analysis.DMS_id == c]
        tmp_d = dict(zip(tmp_df['variable'], tmp_df['value']))
        tmp_d_joint_only = dict(zip(tmp_df[tmp_df['variable'] == 'm12_lrm_pred']['variable'],
                                    tmp_df[tmp_df['variable'] == 'm12_lrm_pred']['value']))
        auc_d[c] = tmp_d
        c_flt_tmp.loc[c_flt_tmp['DMS_id'] == c, 'AUC'] = c_flt_tmp.loc[c_flt_tmp['DMS_id'] == c, 'level_1'].map(
            auc_d[c])
        auc_df = pd.DataFrame.from_dict(auc_d, orient='index')
        auc_df = auc_df.reset_index().melt(value_vars=['m12_lrm_pred'], id_vars=['index'])
        auc_df = auc_df.rename(columns={'index': 'DMS_id', 'variable': 'level_1', 'value': 'AUC'})
        auc_dfs.append(auc_df)
    final_df = pd.concat([c_flt_tmp, pd.concat(auc_dfs)])
    print(final_df.level_1.unique())
    final_df_AUC = final_df.dropna(subset=['AUC']).drop_duplicates().sort_values(by='DMS_id').reset_index()
    print(final_df_AUC.level_1.unique())
    final_df_m12_ols_pred = final_df.dropna(subset=['m12_ols_pred']).drop_duplicates().sort_values(
        by='DMS_id').reset_index()
    print(final_df_m12_ols_pred.level_1.unique())
    final_df = pd.concat([final_df_AUC, final_df_m12_ols_pred])
    final_df['pae'] = final_df['pae'].fillna(method="ffill")
    final_df['MSA_N_eff'] = final_df['MSA_N_eff'].fillna(method="ffill")
    print(final_df.level_1.unique())
    return final_df.drop_duplicates()


# Comparing performance in terms of AUC
def get_auc_df(df_merged, MSA_N_eff_d, pae_d_dms, mean_pLDDT_map):
    df_merged_with_bins = df_merged.dropna(subset=['DMS_score_bin'])

    c1 = 0
    c2 = 0
    c12 = 0
    for id in df_merged_with_bins.DMS_id.unique():
        from sklearn.metrics import roc_auc_score
        auc1 = roc_auc_score(df_merged_with_bins[df_merged_with_bins.DMS_id == id]['DMS_score_bin'].astype(int),
                             df_merged_with_bins[df_merged_with_bins.DMS_id == id]['m1_lrm_pred'])  # evorator_scannet
        auc2 = roc_auc_score(df_merged_with_bins[df_merged_with_bins.DMS_id == id]['DMS_score_bin'].astype(int),
                             df_merged_with_bins[df_merged_with_bins.DMS_id == id]['m2_lrm_pred'])  # tranception eve
        auc12 = roc_auc_score(df_merged_with_bins[df_merged_with_bins.DMS_id == id]['DMS_score_bin'].astype(int),
                              df_merged_with_bins[df_merged_with_bins.DMS_id == id]['m12_lrm_pred'])  # joint
        m = ['evorator_scannet', 'tranception_eve', 'joint']
        aucs = [auc1, auc2, auc12]
        best_auc = max(aucs)
        best_auc_idx = aucs.index(best_auc)
        best_method = m[best_auc_idx]
        print('best', best_method, best_auc)
        if best_method == m[0]:
            c1 += 1
        elif best_method == m[1]:
            c2 += 1
        else:
            c12 += 1

        df_merged_with_bins.loc[df_merged_with_bins.DMS_id == id, 'auc1'] = auc1
        df_merged_with_bins.loc[df_merged_with_bins.DMS_id == id, 'auc2'] = auc2
        df_merged_with_bins.loc[df_merged_with_bins.DMS_id == id, 'auc12'] = auc12

    df_merged_with_bins_for_analysis = pd.melt(df_merged_with_bins, value_vars=['auc1', 'auc2', 'auc12'],
                                               id_vars=['DMS_id'])
    df_merged_with_bins_for_analysis = df_merged_with_bins_for_analysis.drop_duplicates()
    df_merged_with_bins_for_analysis['MSA_N_eff'] = df_merged_with_bins_for_analysis['DMS_id'].map(MSA_N_eff_d)

    df_merged_with_bins_for_analysis['pae'] = df_merged_with_bins_for_analysis['DMS_id'].map(pae_d_dms)
    df_merged_with_bins_for_analysis['mean_pLDDT'] = df_merged_with_bins_for_analysis['DMS_id'].map(mean_pLDDT_map)

    return df_merged_with_bins_for_analysis


def get_cor_m(df_meta, df, MSA_N_eff_d, pae_d_dms, pLDDT_d):
    df_merged = pd.merge(df_meta, df, how='left', left_on='UniProt_ID', right_on='dms_id')

    df_merged = reduce_memory_usage.reduce_memory_usage(df_merged)
    t1 = df_merged.melt(value_vars=['evorator_pred', 'Tranception_L_no_retrieval', 'Tranception_S_retrieval',
                                    'Tranception_M_retrieval', 'Tranception_L_retrieval', 'EVE_single',
                                    'EVE_ensemble', 'MSA_Transformer_single', 'MSA_Transformer_ensemble',
                                    'ESM1v_single', 'ESM1v_ensemble', 'Wavenet', 'DeepSequence_single',
                                    'DeepSequence_ensemble', 'Site_Independent', 'EVmutation', 'RITA_s',
                                    'RITA_m', 'RITA_l', 'RITA_xl', 'RITA_ensemble', 'Progen2_small',
                                    'Progen2_medium', 'Progen2_base', 'Progen2_large', 'Progen2_xlarge',
                                    'Progen2_ensemble', 'Ensemble_Tranception_EVE',
                                    'blosum_pred'], id_vars=['DMS_filename', 'UniProt_ID', 'DMS_score'],
                        value_name='Predicted_score', var_name='Method').dropna()
    import seaborn as sns
    t2 = t1[t1.Method.isin(['Ensemble_Tranception_EVE', 'evorator_pred'])]
    # sns.lmplot(x='Predicted_score',y='DMS_score',data=t2[t2.UniProt_ID==t1.UniProt_ID.values[0]],hue='Method',lowess=True)
    df_merged['targed_seq_len'] = df_merged['DMS_id'].map(
        df_merged.DMS_id.value_counts()[df_merged.DMS_id.value_counts() > 1].to_dict())
    df_merged_grouped = df_merged.groupby('DMS_id')
    try:
              a = df_merged_grouped[['evorator_pred', 'Tranception_L_no_retrieval', 'Tranception_S_retrieval',
                              'Tranception_M_retrieval', 'Tranception_L_retrieval', 'EVE_single',
                              'EVE_ensemble', 'MSA_Transformer_single', 'MSA_Transformer_ensemble',
                              'ESM1v_single', 'ESM1v_ensemble', 'Wavenet', 'DeepSequence_single',
                              'DeepSequence_ensemble', 'Site_Independent', 'EVmutation', 'RITA_s',
                              'RITA_m', 'RITA_l', 'RITA_xl', 'RITA_ensemble', 'Progen2_small',
                              'Progen2_medium', 'Progen2_base', 'Progen2_large', 'Progen2_xlarge',
                              'Progen2_ensemble', 'Ensemble_Tranception_EVE',
                              'blosum_pred', 'm1_ols_pred', 'm2_ols_pred', 'm12_ols_pred', 'DMS_score']].corr(
              method='spearman').abs().dropna()
       except Exception as e:
              print(e)
              a = df_merged_grouped[['evorator_pred', 'Tranception_L_no_retrieval', 'Tranception_S_retrieval',
                                     'Tranception_M_retrieval', 'Tranception_L_retrieval', 'EVE_single',
                                     'EVE_ensemble', 'MSA_Transformer_single', 'MSA_Transformer_ensemble',
                                     'ESM1v_single', 'ESM1v_ensemble', 'Wavenet', 'DeepSequence_single',
                                     'DeepSequence_ensemble', 'Site_Independent', 'EVmutation', 'RITA_s',
                                     'RITA_m', 'RITA_l', 'RITA_xl', 'RITA_ensemble', 'Progen2_small',
                                     'Progen2_medium', 'Progen2_base', 'Progen2_large', 'Progen2_xlarge',
                                     'Progen2_ensemble', 'Ensemble_Tranception_EVE',
                                     'blosum_pred', 'DMS_score']].corr(
                  method='spearman').abs().dropna()
    c = a
    c = c.reset_index()
    import seaborn as sns
    c_flt = c[c.DMS_score > 0]
    c_flt_tmp = c_flt[c_flt.level_1 != 'DMS_score']
    # c_flt_tmp[c_flt_tmp.level_1.isin(['evorator_pred','Ensemble_Tranception_EVE'])]['targret_seq_len']=c_flt_tmp[c_flt_tmp.level_1.isin(['evorator_pred','Ensemble_Tranception_EVE'])]['DMS_id'].map(df_merged.DMS_id.value_counts()[df_merged.DMS_id.value_counts()>1].to_dict())
    # MSA_N_eff_d = dict(zip(df_meta.DMS_id, df_meta.MSA_N_eff))
    c_flt_tmp['MSA_N_eff'] = c_flt_tmp['DMS_id'].map(MSA_N_eff_d)
    c_flt_tmp['mean_pLDDT'] = c_flt_tmp['DMS_id'].map(mean_pLDDT_map)
    df_meta2['pae'] = df_meta2['Entry'].map(pae_d)
    df_meta['pae'] = df_meta['UniProt_ID'].map(pae_d_dms)
    c_flt_tmp['pae'] = c_flt_tmp['DMS_id'].map(pae_d_dms)

    return c_flt_tmp


##
# p0 = "FINAL_DMS_RESULTS_HUANG_EVORATOR_SCANNET.csv" > output for this input to r to get p
p = "FINAL_DMS_RESULTS_EVORATOR_SCANNET_BIVARIATE.csv"
# p2 = "FINAL_DMS_RESULTS_OLD_EVORATOR.csv"   # output for this input to r to get p3
p3 = "FINAL_DMS_RESULTS_EVORATOR_OLD_BIVARIATE.csv"
# p = "individual_vs_bivar_analysis.csv"
df = dt.fread(p).to_pandas()
df_old = dt.fread(p3).to_pandas()
meta_data = "ProteinGym_reference_file_substitutions.csv"
meta_data2 = "proteingym_to_uniprot_ids.tsv"
df_meta = dt.fread(meta_data).to_pandas()
df_meta2 = dt.fread(meta_data2).to_pandas()

MSA_N_eff_d = dict(zip(df_meta.DMS_id, df_meta.MSA_N_eff))
import pickle

with open('uniprot_to_pae_d.pkl', 'rb') as f:
    pae_d = pickle.load(f)
with open('mean_pLDDT_map.data', 'rb') as f:
    mean_pLDDT_map = pickle.load(f)

df_meta2['pae'] = df_meta2['Entry'].map(pae_d)
pae_d_dms = dict(zip(df_meta2['From'], df_meta2['pae']))
df_meta['pae'] = df_meta['UniProt_ID'].map(pae_d_dms)
pae_d_dms = dict(zip(df_meta.DMS_id, df_meta.pae))

df_meta2['mean_pLDDT'] = df_meta2['Entry'].map(mean_pLDDT_map)
mean_pLDDT_map = dict(zip(df_meta2['From'], df_meta2['mean_pLDDT']))
df_meta['mean_pLDDT'] = df_meta['UniProt_ID'].map(mean_pLDDT_map)
mean_pLDDT_map = dict(zip(df_meta.DMS_id, df_meta.mean_pLDDT))

c_flt_tmp = get_cor_m(df_meta, df, MSA_N_eff_d, pae_d_dms, mean_pLDDT_map)
c_flt_tmp_old = get_cor_m(df_meta, df_old, MSA_N_eff_d, pae_d_dms, mean_pLDDT_map)
df_merged = pd.merge(df_meta, df, how='left', left_on='UniProt_ID', right_on='dms_id')
df_merged_old = pd.merge(df_meta, df_old, how='left', left_on='UniProt_ID', right_on='dms_id')
df_merged_with_bins_for_analysis = get_auc_df(df_merged, MSA_N_eff_d, pae_d_dms, mean_pLDDT_map)
df_merged_old_with_bins_for_analysis = get_auc_df(df_merged_old, MSA_N_eff_d, pae_d_dms, mean_pLDDT_map)
##
final_df_new = get_final_df(df_merged_with_bins_for_analysis, c_flt_tmp)
final_df_old = get_final_df(df_merged_old_with_bins_for_analysis, c_flt_tmp_old)
final_df_old = final_df_old.replace('evorator_pred', 'EvoRator')
final_df = pd.concat([final_df_new, final_df_old[final_df_old.level_1 == 'EvoRator']])
final_df = final_df.replace('evorator_pred', 'EvoRator2')
final_df = final_df.replace('m12_lrm_pred', 'EvoRator2 + Ensemble_Tranception_EVE')
final_df = final_df.replace('m12_ols_pred', 'EvoRator2 + Ensemble_Tranception_EVE')
final_df = final_df[final_df.level_1.isin(
    ['EvoRator', 'EvoRator2', 'Ensemble_Tranception_EVE', 'EvoRator2 + Ensemble_Tranception_EVE'])]
final_df_1 = final_df[final_df.level_1 == 'EvoRator2 + Ensemble_Tranception_EVE'].sort_values('DMS_id')
final_df_1 = final_df_1.groupby('DMS_id').first()[['DMS_score', 'AUC']].reset_index()
final_df_1['level_1'] = 'EvoRator2 + Ensemble_Tranception_EVE'
auc_d_joint = dict(zip(final_df_1['level_1'], final_df_1['AUC']))
rho_d_joint = dict(zip(final_df_1['level_1'], final_df_1['DMS_score']))
final_df = pd.concat(
    [final_df[final_df.level_1.isin(['EvoRator', 'EvoRator2', 'Ensemble_Tranception_EVE'])].drop(columns=['index']),
     final_df_1])
# final_df[final_df.level_1=='EvoRator2 + Ensemble_Tranception_EVE'].sort_values('DMS_id')[['DMS_id','DMS_score','AUC']]['DMS_score'].fillna(method='bfill')
# final_df['DMS_score'] = final_df['DMS_score'].fillna(method='bfill')
# final_df['AUC'] = final_df['AUC'].fillna(method='bfill')
final_df = final_df.sort_values('DMS_id')
final_df['pae'] = final_df['pae'].fillna(method="ffill")
final_df['MSA_N_eff'] = final_df['MSA_N_eff'].fillna(method="ffill")

##
final_df.dropna(axis=1).to_csv('FINAL_DMS_AUC_SPEARMAN_RHO.csv', index=False)
##
##
from scipy import stats

stats.spearmanr(final_df[final_df.level_1 == 'EvoRator2']['AUC'], final_df[final_df.level_1 == 'EvoRator2']['pae'])
stats.spearmanr(final_df[final_df.level_1 == 'Ensemble_Tranception_EVE']['AUC'],
                final_df[final_df.level_1 == 'Ensemble_Tranception_EVE']['MSA_N_eff'])
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

## Comparing performance in terms of Spearman correlation coefficient ##

import matplotlib.pyplot as plt

ax = sns.pointplot(x='DMS_score', y='DMS_id', hue='level_1', data=c_flt_tmp[c_flt_tmp.level_1.isin(
    ['evorator_pred', 'Tranception_L_retrieval', 'EVE_single', 'ESM1v_single', 'MSA_Transformer_single'])].sort_values(
    by=['MSA_Transformer_single'], ascending=False), capsize=.4, join=False)
ax = sns.pointplot(x='DMS_score', y='DMS_id', hue='level_1', data=c_flt_tmp[c_flt_tmp.level_1.isin(
    ['evorator_pred', 'Tranception_L_retrieval', 'EVE_single', 'ESM1v_single', 'MSA_Transformer_single'])].sort_values(
    by=['MSA_Transformer_single'], ascending=False), capsize=.4, join=False)
ax = sns.pointplot(x='DMS_score', y='level_1', data=c_flt_tmp[c_flt_tmp.level_1.isin(
    ['evorator_pred', 'Tranception_L_retrieval', 'EVE_single', 'ESM1v_single', 'MSA_Transformer_single'])].sort_values(
    by=['evorator_pred'], ascending=False), capsize=.4, join=False)
ax = sns.pointplot(x='DMS_score', y='level_1', data=c_flt_tmp[
    c_flt_tmp.level_1.isin(['evorator_pred', 'Ensemble_Tranception_EVE', 'm12_ols_pred'])].sort_values(
    by=['Ensemble_Tranception_EVE'], ascending=False), capsize=.4, join=False)
ax= sns.pointplot(x='DMS_score',y='DMS_id',hue='level_1',data=c_flt_tmp[c_flt_tmp.level_1.isin(['blosum_pred','evorator_pred','Ensemble_Tranception_EVE','m12_ols_pred'])].sort_values(by=['Ensemble_Tranception_EVE'],ascending=False),capsize=.4,join=False)
a=c_flt_tmp[c_flt_tmp.level_1=='m12_ols_pred']['DMS_score'].values - c_flt_tmp[c_flt_tmp.level_1=='Ensemble_Tranception_EVE']['DMS_score'].values
s=sms.CompareMeans(sms.DescrStatsW(c_flt_tmp[c_flt_tmp.level_1=='m12_ols_pred']['DMS_score']),sms.DescrStatsW(c_flt_tmp[c_flt_tmp.level_1=='Ensemble_Tranception_EVE']['DMS_score'])) # CI for the mean difference
# handles, labels = plt.gca().get_legend_handles_labels()
# plt.legend(handles,['EvoRator2','Tranception w/ retrieval','EVE','ESM1v','MSA Transformer'],loc='upper left',frameon=False,bbox_to_anchor=(1,1))
# plt.text(0,2.5,' $\mathregular{R^2}$ = 0.56',size=12, fontweight="bold")
# ax.spines['left'].set_position('zero')
ax.spines['right'].set_color('none')
# ax.spines['bottom'].set_position('zero')
ax.spines['top'].set_color('none')
# ax.xaxis.set_label_coords(1.05, -0.025)
# ax.yaxis.set_label_coords(-0.025, -1.05)
plt.xlabel('Spearman ' + r"$\rho$", size=14, fontweight="bold")
plt.ylabel('')

##
# figure 1
## Comparing performance in terms of Spearman
import matplotlib.pyplot as plt
import seaborn as sns

plt.rcParams["font.family"] = "Times New Roman"

fig = plt.figure(figsize=(8, 9))
ax = sns.pointplot(x='DMS_score', y='DMS_id', hue='level_1',
                   data=final_df.sort_values(by=['DMS_score'], ascending=False), capsize=.4, join=False,
                   plot_kws=dict(alpha=0.3), scale=1)
plt.setp(ax.collections, alpha=.3)  # for the markers

handles, labels = plt.gca().get_legend_handles_labels()
order = [0, 1, 2, 3]
plt.legend([handles[idx] for idx in order],
           [['EvoRator2 +\nEnsemble Tranception EVE', 'Ensemble Tranception EVE', 'EvoRator2', 'EvoRator'][idx] for idx
            in order], frameon=True, bbox_to_anchor=(0.475, 0.2))
ax.spines['right'].set_color('none')
# ax.spines['bottom'].set_position('zero')
ax.spines['top'].set_color('none')
# ax.xaxis.set_label_coords(1.05, -0.025)
# ax.yaxis.set_label_coords(-0.025, -1.05)
plt.xlabel('Spearman ' + r"$\rho$", size=12, fontweight="bold")
plt.ylabel('')
plt.tight_layout()

plt.savefig('fig1.tiff', dpi=300)
##
# figure 2
## Comparing performance in terms of AUC

fig = plt.figure(figsize=(8, 9))

ax = sns.pointplot(x='AUC', y='DMS_id', hue='level_1', data=final_df.sort_values(by=['AUC'], ascending=False),
                   join=False, plot_kws=dict(alpha=0.3), scale=1)
plt.setp(ax.collections, alpha=.3)  # for the markers

handles, labels = plt.gca().get_legend_handles_labels()
order = [0, 1, 2, 3]
plt.legend([handles[idx] for idx in order],
           [['EvoRator2 +\nEnsemble Tranception EVE', 'Ensemble Tranception EVE', 'EvoRator2', 'EvoRator'][idx] for idx
            in order], frameon=True, bbox_to_anchor=(0.575, 0.15))
ax.spines['right'].set_color('none')
# ax.spines['bottom'].set_position('zero')
ax.spines['top'].set_color('none')
# ax.xaxis.set_label_coords(1.05, -0.025)
# ax.yaxis.set_label_coords(-0.025, -1.05)
plt.xlabel('AUC', size=12, fontweight="bold")
plt.ylabel('')
plt.tight_layout()
plt.savefig('fig2.tiff', dpi=300)
##
df_merged_with_bins_for_analysis.to_csv('AUC_vs_MSA_pae.csv',index=False)
c_flt_tmp[c_flt_tmp.level_1.isin(['evorator_pred','Ensemble_Tranception_EVE','m12_ols_pred'])].to_csv('Spearman_vs_MSA_pae.csv',index=False)
##



