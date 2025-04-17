
from sklearn.cluster import KMeans
import blosum as bl
import numpy as np
import datatable as dt
import  seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix, classification_report
from scipy import stats
from statsmodels.regression import linear_model

def PWM(pred,mut):
    aa_freqs = {'A':8.25/100, 'Q':3.93/100, 'L':9.65/100, 'S':6.63/100, 'R':5.53/100,
                'E':6.72/100, 'K':5.80/100, 'T':5.35/100,'N':4.05/100,
                'G':7.08/100, 'M':2.41/100,'W':1.09/100, 'D':5.46/100, 'H':2.27/100,
                'F':3.86/100, 'Y':2.92/100, 'C':1.38/100, 'I':5.91/100, 'P':4.73/100, 'V':6.86/100}
    return np.log2((pred + 1e-5) / aa_freqs[mut])
    # return pred  / aa_freqs[mut]

def analyze_performance(df):
    mat = bl.BLOSUM(62)
    aa = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    aa_pred = [AA + '_pred' for AA in aa]
    aa_pwm_predicted= [AA + "_pwm_predicted" for AA in aa]
    aa_pwm_observed = [AA + "_pwm_observed" for AA in aa]
    aa_pred_baseline = [AA + '_blosum' for AA in aa_pred]

    aa_pred_plot = [AA + '*' for AA in aa]
    if isinstance(df,str):
        df= dt.fread(df).to_pandas()

    print(df.columns.tolist())
    print(df[aa_pred].values)
    print(df[aa].values)

    print(stats.spearmanr(df[aa_pred].values.ravel(),df[aa].values.ravel(),nan_policy='omit'))
    print(stats.pearsonr(df[aa_pred].values.ravel(),df[aa].values.ravel()))
    aa_freq_d = df.pdb_aa.value_counts(normalize=True).to_dict()
    print(aa_freq_d)
    aa_freq_k = aa_freq_d.keys()
    aa_freq_v = list(aa_freq_d.values())
    n_aa = sum(aa_freq_v)
    # print(df[aa_pred])
    # print(df[aa_pred].idxmax(axis=1).str[0])
    df['aa_with_max_predicted_freq'] = df[aa_pred].idxmax(axis=1).str[0]

    for i in range(len(aa_pred_baseline)):
        df[aa_pred_baseline[i]] = df.pdb_aa+aa_pred_baseline[i][0]
        df[aa_pred_baseline[i]] = df[aa_pred_baseline[i]].map(mat)

    for mut in aa_pred:
        df[mut[0] + "_pwm_predicted"] = df[mut].apply(PWM, mut=mut[0])
    for mut in aa_pred:
        df[mut[0] + "_pwm_observed"] = df[mut[0]].apply(PWM, mut=mut[0])

    print(stats.spearmanr(df[aa_pwm_observed].values.ravel(), df[aa_pwm_predicted].values.ravel()))
    print(stats.spearmanr(df[aa_pwm_observed].values.ravel(), df[aa_pred_baseline].values.ravel()))
    print(stats.pearsonr(df[aa_pwm_observed].values.ravel(), df[aa_pwm_predicted].values.ravel()))
    print(stats.pearsonr(df[aa_pwm_observed].values.ravel(), df[aa_pred_baseline].values.ravel()))
    ''' Does it predict conservation well? '''
    ''' How different is it from a constant predictor (I.e. a BLOSUM substitution matrix)'''
    ''' Accuracy vs sequence length '''
    ''' Accuracy vs baseline frequency '''
    ''' How well it predicts the most frequent substitution '''
    '''- Can it predict « unexpected » substitutions (say, hydrophobic to charged residue?). Generally speaking, do the predictions span across all possible PSSMs, or does it make « safe » predictions only (like hydrophobic to hydrophobic)'''
    print(np.diag(df[[*aa_pwm_observed, *aa_pwm_predicted]].corr(method='spearman').iloc[:20,20:]))
    print(np.diag(df[[*aa_pwm_observed, *aa_pred_baseline]].corr(method='spearman').iloc[:20,20:]))
    # plt.show()
    # plt.matshow(df[[*aa_pwm_observed,*aa_pred_baseline]].corr(method='pearson').iloc[:20,20:])
    # fig, ax = plt.subplots(nrows=1, ncols=2,figsize=(6, 4))
    # hb = ax[0].hexbin(df[aa_pwm_predicted].values.ravel(), df[aa_pwm_observed].values.ravel(),
    #                gridsize=20, edgecolors='grey',
    #                cmap='inferno', mincnt=1)
    # ax[0].set_axisbelow(True)
    # ax[0].set_xlabel('Predicted conservation score')
    # ax[0].set_ylabel('Observed conservation score')
    # cb = fig.colorbar(hb, ax=ax[0])
    # ax[0].plot([0, 1], [0, 1], transform=ax[0].transAxes)
    #
    # hb = ax[1].hexbin(df[aa_pred_baseline].values.ravel(), df[aa_pwm_observed].values.ravel(),
    #                gridsize=20, edgecolors='grey',
    #                cmap='inferno', mincnt=1)
    # ax[1].set_axisbelow(True)
    # ax[1].set_xlabel('Predicted conservation score')
    # ax[1].set_ylabel('Observed conservation score')
    # cb = fig.colorbar(hb, ax=ax[1])
    # ax[1].plot([0, 1], [0, 1], transform=ax[1].transAxes)
    # plt.show()
    #
    # plt.matshow(df[[*aa_pwm_observed,*aa_pwm_predicted]].corr(method='pearson').iloc[:20,20:])
    # plt.show()
    # plt.matshow(df[[*aa_pwm_observed,*aa_pred_baseline]].corr(method='pearson').iloc[:20,20:])
    # plt.show()
    # plt.matshow(df[[*aa,*aa_pred]].corr(method='pearson').iloc[:20,20:])
    # plt.show()
    df[aa] = df[aa] + 0.001
    g = df[aa_pred]
    df['aa_with_second_largest_predicted_freq'] = g.mask(g.eq(g.max(axis=1), axis=0) & g.apply(lambda x: ~x.duplicated(), axis=1)).idxmax(axis=1).str[0]

    df['aa_with_max_observed_freq'] = df[aa].idxmax(axis=1).str[0]
    g = df[aa]
    df['aa_with_second_largest_observed_freq'] = g.mask(g.eq(g.max(axis=1), axis=0) & g.apply(lambda x: ~x.duplicated(), axis=1)).idxmax(axis=1).str[0]

    df['aa_with_max_predicted_freq_base'] = df[aa_pred_baseline].idxmax(axis=1).str[0]
    g = df[aa_pred_baseline]
    df['aa_with_second_largest_predicted_freq_base'] = g.mask(g.eq(g.max(axis=1), axis=0) & g.apply(lambda x: ~x.duplicated(), axis=1)).idxmax(axis=1).str[0]

    # print(df['aa_with_max_predicted_freq'])
    # print(df['aa_with_second_largest_predicted_freq'])
    # print(df['aa_with_max_observed_freq'])
    # print(df['aa_with_second_largest_observed_freq'])
    # exit()
    df['match_original_aa_bool'] = df.pdb_aa==df['aa_with_max_predicted_freq'] #  the aa in the primary sequence matches the aa with the largest predicted frequency
    df['match_max_observed_freq_aa_bool'] = df.aa_with_max_observed_freq==df['aa_with_max_predicted_freq'] #  the aa with the largest observed frequency matches the aa with the largest predicted frequency
    df['match_second_largest_observed_freq_aa_bool'] = df.aa_with_second_largest_observed_freq==df['aa_with_second_largest_predicted_freq'] #  the aa with the second largest observed frequency matches the aa with the second largest predicted frequency

    print(classification_report(df.pdb_aa.values,df['aa_with_max_observed_freq'].values,labels=aa))
    print(classification_report(df.pdb_aa.values,df['aa_with_max_predicted_freq'].values,labels=aa))
    print(classification_report(df.pdb_aa.values,df['aa_with_max_predicted_freq_base'].values,labels=aa))
    print(classification_report(df.aa_with_max_observed_freq.values,df['aa_with_max_predicted_freq'].values,labels=aa))
    print(classification_report(df.aa_with_max_observed_freq.values,df['aa_with_max_predicted_freq_base'].values,labels=aa))
    print(classification_report(df.aa_with_second_largest_observed_freq.values,df['aa_with_second_largest_predicted_freq'].values,labels=aa))
    print(classification_report(df.aa_with_second_largest_observed_freq.values,df['aa_with_second_largest_predicted_freq_base'].values,labels=aa))

    df_grouped =  df.groupby('split_by_CATH')[['predicted_pssm_conservation_score','pssm_conservation_score']].corr(method='spearman').reset_index()
    from sklearn.linear_model import LinearRegression
    lr = LinearRegression()
    t=df[['predicted_pssm_conservation_score','pssm_conservation_score','aa_pred_baseline']].fillna(df['pssm_conservation_score'].mean())
    u, s = t['predicted_pssm_conservation_score'].mean(), t['predicted_pssm_conservation_score'].std()
    t['predicted_pssm_conservation_score_std'] = t['predicted_pssm_conservation_score'].apply(lambda x: (x-u)/s)
    lr.fit(t[['predicted_pssm_conservation_score_std',*aa_pred_baseline]],t[['pssm_conservation_score']])

    print(lr.coef_)
    print(lr.intercept_)
    print(lr.score(t[['predicted_pssm_conservation_score_std']],t[['pssm_conservation_score']]))

    lr = LinearRegression()
    t=df[['evorator_conservation_score','pssm_conservation_score']].fillna(df['pssm_conservation_score'].mean())
    u, s = t['evorator_conservation_score'].mean(), t['evorator_conservation_score'].std()
    t['evorator_conservation_score_std'] = t['evorator_conservation_score'].apply(lambda x: (x-u)/s)
    lr.fit(t[['evorator_conservation_score_std']],t[['pssm_conservation_score']])

    print(lr.coef_)
    print(lr.intercept_)
    print(lr.score(t[['evorator_conservation_score_std']],t[['pssm_conservation_score']]))

    exit()
    df_grouped = df_grouped[df_grouped['level_1']=='predicted_pssm_conservation_score']
    sns.violinplot(df_grouped['pssm_conservation_score'],cut=0)
    plt.show()

    for i in aa:
        aa_freq_v[aa.index(i)] = aa_freq_d[i]/n_aa
    print(aa_freq_v)
    plt.matshow(df[[*aa_pred,*aa]].corr(method='spearman').iloc[:20,20:])
    # g = sns.clustermap(df[[*aa_pred,*aa]].corr(method='spearman').iloc[:20,20:])
    print(np.diag(df[[*aa_pred,*aa]].corr(method='spearman').iloc[:20,20:]))
    plt.xticks(list(range(20)),labels=aa)
    plt.yticks(list(range(20)),labels=aa_pred_plot)
    plt.show()
    fig, ax = plt.subplots(figsize=(6, 4))
    hb = ax.hexbin(df['predicted_pssm_conservation_score'], df['pssm_conservation_score'],
                   gridsize=20, edgecolors='grey',
                   cmap='inferno', mincnt=1)
    ax.set_axisbelow(True)
    ax.set_xlabel('Predicted conservation score')
    ax.set_ylabel('Observed conservation score')
    cb = fig.colorbar(hb, ax=ax)
    ax.plot([0, 1], [0, 1], transform=ax.transAxes)
    plt.show()



    sns.jointplot(x=aa_freq_v, y=np.diag(df[[*aa_pred,*aa]].corr(method='spearman').iloc[:20,20:]), kind="reg",color="#4CB391")
    plt.show()
    for i in range(len(aa)):
        sns.regplot(x=df[aa_pred[i]],y=df[aa[i]],data=df, lowess=True, scatter=True)
        plt.show()


print('loading data')
# df= dt.fread(r"results\model_selection\deep_learning_for_proteinnet\final_output_deep_learning_for_proteinnet_n_hidden_30_test_pred.csv").to_pandas()
# df= dt.fread(r"debug_milion_deep_learning_for_proteinnet_n_hidden_30_test_pred.csv").to_pandas()
df= dt.fread(r"C:\Users\natan\Documents\EvolutionPrediction\results\model_selection\deep_learning_for_proteinnet\pp_deep_learning_for_proteinnet_n_hidden_30_test_pred.csv").to_pandas()
print('data loaded')

print(df.shape)
# df = df.iloc[:,:5000]
# df= pd.read_csv(r"results\model_selection\deep_learning_for_proteinnet\final_output_deep_learning_for_proteinnet_n_hidden_30_test_pred.csv")
analyze_performance(df)