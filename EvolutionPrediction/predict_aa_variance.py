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
import click


import itertools
import time
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestRegressor
from sklearn.neural_network import MLPRegressor
from sklearn.svm import SVR
from sklearn.linear_model import LinearRegression, Lasso, Ridge, LassoCV
import pandas as pd
import random
from sklearn.pipeline import make_pipeline
from sklearn.model_selection import cross_val_score, cross_val_predict, train_test_split, RepeatedKFold, GroupKFold, \
    KFold, GridSearchCV, GroupShuffleSplit, permutation_test_score
from sklearn.impute import SimpleImputer
import os
from sklearn.metrics import r2_score
import logging
import pickle
from sklearn.preprocessing import LabelEncoder, OneHotEncoder, StandardScaler
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import make_pipeline, Pipeline
from sklearn.ensemble import RandomForestRegressor, ExtraTreesRegressor
from sklearn.feature_selection import RFECV, VarianceThreshold, RFE, SelectFromModel
import scipy
from sklearn import clone
import training_constants
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegressionCV, LogisticRegression
from sklearn.metrics import log_loss, roc_auc_score
from sklearn.multiclass import OneVsRestClassifier
from sklearn.naive_bayes import MultinomialNB
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import make_pipeline, Pipeline
from sklearn.model_selection import cross_val_score, KFold, cross_val_predict, train_test_split
from sklearn.preprocessing import LabelEncoder, OneHotEncoder, StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import MultiLabelBinarizer
from sklearn.neural_network import MLPClassifier

from tensorflow.python.keras.layers import Input, Dense, Activation, Flatten, Conv2D, BatchNormalization
from tensorflow.python.keras.wrappers.scikit_learn import KerasRegressor
from tensorflow.python.keras.models import Model, load_model
from tensorflow.python.keras.regularizers import l2, l1
from Bio.SubsMat import MatrixInfo as matlist
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.multioutput import MultiOutputClassifier
from tensorflow.python.keras.wrappers.scikit_learn import  KerasClassifier
import pickle
# import lightgbm as lgb
import reduce_memory_usage

data = r"D:\university\projects\nathan\evorator\evorator4_aa_variety_pred.csv"

def NNModel_class(input_dim=None, num_classes=2, n_hidden=20,w=None):
    X_input = Input(input_dim)
    X = Dense(input_dim[0], activation='relu', name='fc2',kernel_regularizer=l2(0.0005))(X_input)
    X = BatchNormalization(name = 'bn0')(X)
    X = Dense(input_dim[0], activation='relu', name='fc3',kernel_regularizer=l2(0.0005))(X)
    X = BatchNormalization(name = 'bn1')(X)
    X = Dense(num_classes, activation='softmax', name='fc_out')(X) # as in preprint
    model = Model(inputs=X_input, outputs=X)
    model.compile('adam', 'kullback_leibler_divergence', metrics=['kullback_leibler_divergence'])  # as in preprint

    return model

@click.command()
@click.argument('test_data', type=str)
@click.argument('model_path', type=str)
@click.argument('preprocessor', type=str)
@click.argument('training_dim', type=int) # training_dim: 467 - 1e5 data, 515 - all data
@click.option('--output-name','-o', type=str, default='final_output')
@click.option('--adress-db','-db', is_flag=True)
def map_residue_variety_to_features(test_data,model_path,preprocessor,training_dim,output_name,adress_db):

    # TEST_DATA_PATH = r"C:\Users\natan\Documents\missense_pathogenecity\source\PDB1BE9_features_and_predictions.csv"
    # TEST_DATA_PATH = r"C:\Users\natan\Documents\missense_pathogenecity\data\1AXB_features_and_predictions.csv"
    # TEST_DATA_PATH = r"C:\Users\natan\Documents\missense_pathogenecity\data\PDB4ZMJ_features_and_predictions.csv"
    # TEST_DATA_PATH = r"C:\Users\natan\Documents\missense_pathogenecity\data\1GNX_features_and_predictions.csv"
    # TEST_DATA_PATH = r"C:\Users\natan\Documents\deg_project\1D66_features_and_predictions.csv"
    # TEST_DATA_PATH = r"C:\Users\natan\Documents\missense_pathogenecity\data\4LVH_features_and_predictions.csv"
    # TEST_DATA_PATH = r"C:\Users\natan\Documents\missense_pathogenecity\data\1AH6_features_and_predictions.csv"
    # TEST_DATA_PATH = r"C:\Users\natan\Documents\missense_pathogenecity\data\1ND4_features_and_predictions.csv"
    # TEST_DATA_PATH = r"C:\Users\natan\Documents\missense_pathogenecity\data\3UBT_features_and_predictions.csv"
    # TEST_DATA_PATH = r"C:\Users\natan\Documents\missense_pathogenecity\data\5M3H_features_and_predictions.csv"
    # TEST_DATA_PATH = r"C:\Users\natan\Documents\missense_pathogenecity\data\6R5K_features_and_predictions.csv"
    # TEST_DATA_PATH = r"C:\Users\natan\Documents\missense_pathogenecity\data\AF-Q9ES00-F1-MODEL_V2-UBE4B-MOUSE_features_and_predictions.csv"
    # TEST_DATA_PATH = r"C:\Users\natan\Documents\missense_pathogenecity\data\4REX_features_and_predictions.csv"
    # TEST_DATA_PATH = r"C:\Users\natan\Documents\missense_pathogenecity\data\2UXY_features_and_predictions.csv"
    # TEST_DATA_PATH = r"C:\Users\natan\Documents\missense_pathogenecity\data\3S4Y_features_and_predictions.csv"
    # TEST_DATA_PATH = r"C:\Users\natan\Documents\missense_pathogenecity\data\2GRN_features_and_predictions.csv"
    # TEST_DATA_PATH = r"C:\Users\natan\Documents\missense_pathogenecity\data\1WYW_features_and_predictions.csv"
    # TEST_DATA_PATH = r"C:\Users\natan\Documents\missense_pathogenecity\data\1CLL_features_and_predictions.csv"
    # TEST_DATA_PATH = r"C:\Users\natan\Documents\missense_pathogenecity\data\6EZM_features_and_predictions.csv"
    # TEST_DATA_PATH = r"C:\Users\natan\Documents\missense_pathogenecity\data\AF-Q9ES00-F1-model_v2-IF1-ecoli_features_and_predictions.csv"
    # TEST_DATA_PATH = r"C:\Users\natan\Documents\missense_pathogenecity\data\121P_features_and_predictions.csv"
    # TEST_DATA_PATH = r"C:\Users\natan\Documents\missense_pathogenecity\data\AF-Q9ES00-F1-MODEL_V2-RL401-YEAST_features_and_predictions.csv"
    # TEST_DATA_PATH = r"C:\Users\natan\Documents\missense_pathogenecity\data\1PME_features_and_predictions.csv"
    # TEST_DATA_PATH = r"C:\Users\natan\Documents\missense_pathogenecity\data\4YH5_features_and_predictions.csv"
    # TEST_DATA_PATH = r"C:\Users\natan\Documents\missense_pathogenecity\data\1D5R_features_and_predictions.csv"
    # TEST_DATA_PATH = r"C:\Users\natan\Documents\missense_pathogenecity\data\2BZG_features_and_predictions.csv"
    # TEST_DATA_PATH = r"C:\Users\natan\Documents\missense_pathogenecity\data\1A53_features_and_predictions.csv"
    # TEST_DATA_PATH = r"C:\Users\natan\Documents\missense_pathogenecity\data\1I4N_features_and_predictions.csv"
    # TEST_DATA_PATH = r"C:\Users\natan\Documents\missense_pathogenecity\data\1VC4_features_and_predictions.csv"
    # TEST_DATA_PATH = r"C:\Users\natan\Documents\missense_pathogenecity\data\7LYB_features_and_predictions.csv"
    # TEST_DATA_PATH = r"C:\Users\natan\Documents\missense_pathogenecity\data\1JNX_features_and_predictions.csv"
    # TEST_DATA_PATH = r"C:\Users\natan\Documents\conservation_paper\revision_submission\final_submission_revision\2YMK_features_and_predictions.csv"
    # TEST_DATA_PATH = r"C:\Users\natan\Documents\missense_pathogenecity\output\evorator4_aa_variety_pred.csv"
    #


    try:

        print('ATTEMPTING TO LOAD THE MODEL USING LOAD MODEL')
        # if laptop:

        model = load_model(model_path)
    except Exception as e:
            print('FAILED TO LOAD MODEL')
            print(e)
            try:
                print('TRYING A DIFFERENT APPROACH TO LOAD MODEL')
                model = NNModel_class(input_dim=[training_dim], num_classes=20,
                                      n_hidden=147)  #
                model.load_weights(model_path)
            except Exception as e:
                print('FAILED TO LOAD MODEL')
                print(e)
                return

        # except:
        #     TRAINED_MODEL = saved_model.load(MODEL_PATH)
    print('MODEL LOADED')
    if not laptop:
        OUTDIR = r'/bioseq/evorator/EvolutionPrediction/substitutions/evorator_scannet_predictions'
    else:
        OUTDIR = r'substitutions/evorator_scannet_predictions'
    PRED_OUTPATH = os.path.join(OUTDIR,output_name+'.csv')
    # TRAINED_MODEL = pickle.load(open(model))
    import sklearn
    # print(sklearn.__version__)
    PREPROCESSOR = pickle.load(open(preprocessor,'rb'))
    matrix = matlist.blosum62
    test_data = pd.read_csv(test_data)
    try:
        ppi_map = pickle.load(open(r"/bioseq/evorator/consurf_db_random/ppi_map.data",'rb'))
        test_data['protein_interaction'] = test_data['pdb_position'].map(ppi_map)
    except Exception as e:
        print(e)

    test_data = reduce_memory_usage.reduce_memory_usage(test_data)
    aa = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    # print(f'test data dimensions = {test_data.shape}')
    test_data = test_data.rename(columns={'total_ss_loop_neigh': 'total_ss_nan_neigh'})
    # print(test_data['pdb_aa'].isna().sum())
    if adress_db:
        test_data['pdb_aa'] = test_data['pdb_aa'].fillna(test_data['pdb_aa_x'].fillna(test_data['pdb_aa_y']))
        test_data = test_data.dropna(subset=aa)  # complete case analysis

    test_data = test_data.rename(columns={'total_ss_loop_neigh': 'total_ss_nan_neigh'})
    test_data = test_data.rename(columns={'predicted_score': 'predicted_evorator_score'})
    blosum_d = {}
    for AA1 in aa:
        tmp_vec = []
        for AA2 in aa:
            tmp_vec.append(matrix.get((AA1, AA2), matrix.get((AA2, AA1))))
        #     # print(AA2,sep=' ')
        # print(len(tmp_vec))
        blosum_d[AA1] = tmp_vec
    # print(blosum_d)
    # print(test_data['pdb_aa'].astype(str))
    # print(test_data['pdb_aa'].astype(str).map(blosum_d))
    test_data[['blosum'+str(i) for i in range(len(blosum_d['A']))]] = test_data['pdb_aa'].map(blosum_d).apply(pd.Series)
    aa_freqs = {'A': 8.25 / 100, 'Q': 3.93 / 100, 'L': 9.65 / 100, 'S': 6.63 / 100, 'R': 5.53 / 100,
                'E': 6.72 / 100, 'K': 5.80 / 100, 'T': 5.35 / 100, 'N': 4.05 / 100,
                'G': 7.08 / 100, 'M': 2.41 / 100, 'W': 1.09 / 100, 'D': 5.46 / 100, 'H': 2.27 / 100,
                'F': 3.86 / 100, 'Y': 2.92 / 100, 'C': 1.38 / 100, 'I': 5.91 / 100, 'P': 4.73 / 100, 'V': 6.86 / 100}
    test_data['aa_freq_ExPASy'] = test_data['pdb_aa'].map(aa_freqs)


    # print(test_data.columns.tolist())

    # not_features = ['index', 'pdb_position', 'pdb_id', 'chain', 'pdb', 'pdb_x', 'pdb_y', 'chain_x', 'chain_y',
    #                 'pdb_aa_y', 'Functional', 'd', 'g','_merge', 'blosum0', 'blosum1', 'blosum2',
    #                 'Unnamed: 0','blosum3', 'blosum4', 'blosum5', 'blosum6', 'blosum7', 'blosum8', 'blosum9', 'blosum10','predicted_evorator_score', 'ResInd', 'Score'
    #                    , 'blosum11', 'blosum12', 'blosum13', 'blosum14', 'blosum15', 'blosum16', 'blosum17', 'blosum18', 'blosum19', 'aa_freq_ExPASy',
    #                 '1/d', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M',
    #                 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', 'pdb_aa_x', 'normalized_score',
    #                 'total_nan_neigh', 'zr4s_JC'] \
    #                + [f for f in test_data.columns.tolist()
    #                   if f.startswith('neighbor_pos_') or f.startswith("n2v")]

    not_features = ['index', 'pdb_position', 'pdb_id', 'chain', 'pdb', 'pdb_x', 'pdb_y', 'chain_x', 'chain_y','Unnamed: 0','ResInd','Score','_merge'
                    'pdb_aa_y', 'Functional', 'd', 'g','blosum3', 'blosum4', 'blosum5', 'blosum6', 'blosum7', 'blosum8', 'blosum9', 'blosum10',
                    '1/d', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M',
                    'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', 'pdb_aa_x', 'normalized_score',
                    'total_nan_neigh', 'zr4s_JC'] \
                   + [f for f in test_data.columns.tolist()
                      if f.startswith('neighbor_pos_') or f.startswith("n2v")]

    categorical_features = ['glycosylation', 'protein_interaction', 'disorder', 'binding', 'catalysis', 'pdb_aa',
                            'aa_group_5', 'aa_group_HP', 'structure']
    '''numeric_transformer_rf = Pipeline(
        steps=[('imputer', SimpleImputer(strategy='median')), ('scaler', StandardScaler()),
               ('variance_filter', VarianceThreshold())])
    categorical_transformer = Pipeline(
        steps=[('imputer', SimpleImputer(strategy='most_frequent')), ('onehot', ohe),
               ('variance_filter', VarianceThreshold())])
    column_transformer = ColumnTransformer(
        transformers=[
            ('num', numeric_transformer_rf, numeric_features),
            ('cat', categorical_transformer, categorical_features),
        ])
    steps = [('preprocessor', column_transformer)]
    # if not os.path.exists(FITTED_PREPROCESSOR_OUTPATH):

    preprocessor = Pipeline(steps)'''

#     # print(PREPROCESSOR.named_steps['preprocessor'].transformers_[1][1]['onehot'].get_feature_names(categorical_features))
#     # print(PREPROCESSOR.named_steps['preprocessor'].transformers_[1][1]['onehot'].transform(test_data[[*categorical_features]]))
    # exit()
#     # print(PREPROCESSOR.named_steps['preprocessor'].transformers_[1][1]['onehot'].get_feature_names(numerical_features))
    all_features = [f for f in test_data.columns.tolist() if f not in not_features]
    # numeric_features = [f for f in all_features if f not in [*categorical_features,'pdb_aa_y']]
    numeric_features = training_constants.numeric_features
    categorical_features = training_constants.categorical_features
    print(f'numeric features={numeric_features}')
    print(f'categorical features={categorical_features}')

    test_data = test_data.rename(columns={'predicted_evorator_score':'evorator_conservation_score',
                                          'total_ss_nan_neigh':'total_ss_loop_neigh'})
    for f in numeric_features:
        if f not in test_data.columns.tolist():
            test_data[f] = np.nan

    X_test = test_data[[*numeric_features, *categorical_features]]  # for comparing with Debora Marks
    # print(test_data.columns.tolist())
    if adress_db:
        y_test = test_data[aa].astype(np.int8)

        # print(f'X_test={X_test.shape},y_test={y_test.shape}')

    # print(f'X_test={X_test.shape}')
    del_v = ['2_clique', '3_clique', '4_clique', '5_clique', '6_clique', '7_clique', 'average_neighbor_degree',
     'betweenness_centrality', 'clustering_coefficient', 'degree_centrality', 'eigenvector_centrality', 'graphlet1',
     'graphlet10', 'graphlet11', 'graphlet12', 'graphlet13', 'graphlet14', 'graphlet15', 'graphlet16', 'graphlet17',
     'graphlet18', 'graphlet19', 'graphlet2', 'graphlet20', 'graphlet21', 'graphlet22', 'graphlet23', 'graphlet24',
     'graphlet25', 'graphlet26', 'graphlet27', 'graphlet28', 'graphlet29', 'graphlet3', 'graphlet30', 'graphlet31',
     'graphlet32', 'graphlet33', 'graphlet34', 'graphlet35', 'graphlet36', 'graphlet37', 'graphlet38', 'graphlet39',
     'graphlet4', 'graphlet40', 'graphlet41', 'graphlet42', 'graphlet43', 'graphlet44', 'graphlet45', 'graphlet46',
     'graphlet47', 'graphlet48', 'graphlet49', 'graphlet5', 'graphlet50', 'graphlet51', 'graphlet52', 'graphlet53',
     'graphlet54', 'graphlet55', 'graphlet56', 'graphlet57', 'graphlet58', 'graphlet59', 'graphlet6', 'graphlet60',
     'graphlet61', 'graphlet62', 'graphlet63', 'graphlet64', 'graphlet65', 'graphlet66', 'graphlet67', 'graphlet68',
     'graphlet69', 'graphlet7', 'graphlet70', 'graphlet71', 'graphlet72', 'graphlet73', 'graphlet8', 'graphlet9',
     'median_rsa_neigh', 'median_wcn_ca_neigh', 'median_wcn_sc_neigh', 'node_degree', 'rsa', 'total_A_neigh',
     'total_Aliphatic_neigh', 'total_Aromatic_neigh', 'total_C_neigh', 'total_Charged_neigh', 'total_D_neigh',
     'total_Diverse_neigh', 'total_E_neigh', 'total_F_neigh', 'total_G_neigh', 'total_H_neigh',
     'total_Hydrophobic_neigh', 'total_I_neigh', 'total_K_neigh', 'total_L_neigh', 'total_M_neigh', 'total_N_neigh',
     'total_P_neigh', 'total_Polar_neigh', 'total_Q_neigh', 'total_R_neigh', 'total_S_neigh', 'total_T_neigh',
     'total_Tiny_neigh', 'total_V_neigh', 'total_W_neigh', 'total_Y_neigh', 'total_catalytic_neigh',
     'total_contact_neigh', 'total_disordered_neigh', 'total_glycosylated_neigh', 'total_interface_neigh',
     'total_non_interfacing_neigh', 'total_site_neigh', 'total_ss_B_neigh', 'total_ss_E_neigh', 'total_ss_G_neigh',
     'total_ss_H_neigh', 'total_ss_I_neigh', 'total_ss_P_neigh', 'total_ss_S_neigh', 'total_ss_T_neigh',
     'total_ss_nan_neigh', 'wcn_ca', 'wcn_sc', 'predicted_evorator_score', 'blosum0', 'blosum1', 'blosum2', 'blosum3',
     'blosum4', 'blosum5', 'blosum6', 'blosum7', 'blosum8', 'blosum9', 'blosum10', 'blosum11', 'blosum12', 'blosum13',
     'blosum14', 'blosum15', 'blosum16', 'blosum17', 'blosum18', 'blosum19', 'aa_freq_ExPASy']
    # if X_test.shape[1] != 164:
#     #     print(f'X_test shape is {X_test.shape[1]}. Expected 164')
#     #     print(f'Dropping total_non_interfacing_neigh')
    #     X_test = X_test.drop(['total_non_interfacing_neigh'],axis=1)

    print('TRANSFORMING TEST FEATURES')
    X_test_proc = PREPROCESSOR.transform(X_test)
    print(f'X test shape {X_test_proc.shape}')
    print('TRYING TO PREDICT')
    print(model)
    print(model.summary())
    pred = model.predict(X_test_proc)
    print('prediction done')

    #@$#@%#@%$#$#%#$#$^$#^#$^$#
    #@$#@%#@%$#$#%#$#$^$#^#$^$#
    #@$#@%#@%$#$#%#$#$^$#^#$^$#
    #@$#@%#@%$#$#%#$#$^$#^#$^$#
    #@$#@%#@%$#$#%#$#$^$#^#$^$#
    #!#!#!#!#!# save fitted preprocessing pipeline to predict new cases #!#!#!#!#!#!#!#!#!#!#!
    #@$#@%#@%$#$#%#$#$^$#^#$^$#
    #@$#@%#@%$#$#%#$#$^$#^#$^$#
    #@$#@%#@%$#$#%#$#$^$#^#$^$#
    #@$#@%#@%$#$#%#$#$^$#^#$^$#

#     # print(len(pred))
#     # print(pred[0].shape)
    # pred_prob = [i[:, 1].tolist() for i in pred]
    if adress_db:
        score = roc_auc_score(y_test,pd.DataFrame(pred_prob).T.values,average=None)
        # print(f'AUC = {score}')
        score = roc_auc_score(y_test,pd.DataFrame(pred_prob).T.values,average='micro')
        # print(f'AUC = {score}')
        # print(np.mean(score))



    # pred_df = pd.DataFrame(pd.DataFrame(pred).T.values,
    #                        columns=['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T',
    #                                 'V', 'W', 'Y'], index=test_data.pdb_position.tolist())
    pred_df = pd.DataFrame(pd.DataFrame(pred).values,
                           columns=aa, index=test_data.pdb_position.tolist())
    pred_df = pd.concat([test_data['pdb_aa'].reset_index(), pred_df.reset_index()], axis=1)
    # stats.spearmanr(pred_df.values.ravel(), pred_training_set.ravel())
    for c in pred_df.columns.tolist():
        if c in aa:
            pred_df.loc[pred_df.pdb_aa==c,'wt_val'] = pred_df[pred_df.pdb_aa==c][c]

    for c in pred_df.columns.tolist():
        if c in aa:
            pred_df[c] = np.log2(pred_df[c]/pred_df['wt_val'])
    print(f'final output = {pred_df}')
    pred_df['index'] = test_data.pdb_position.tolist()
    print(pred_df)
    pred_df.to_csv(PRED_OUTPATH, index=False)
    return PRED_OUTPATH



if __name__ == '__main__':
    import itertools
    import time
    import numpy as np
    import seaborn as sns
    import matplotlib.pyplot as plt
    from sklearn.decomposition import PCA
    from sklearn.ensemble import RandomForestRegressor
    from sklearn.neural_network import MLPRegressor
    from sklearn.svm import SVR
    from sklearn.linear_model import LinearRegression, Lasso, Ridge, LassoCV
    import pandas as pd
    import random
    from sklearn.pipeline import make_pipeline
    from sklearn.model_selection import cross_val_score, cross_val_predict, train_test_split, RepeatedKFold, GroupKFold, \
        KFold, GridSearchCV, GroupShuffleSplit, permutation_test_score
    from sklearn.impute import SimpleImputer
    import os
    from sklearn.metrics import r2_score
    import logging
    import pickle
    from sklearn.preprocessing import LabelEncoder, OneHotEncoder, StandardScaler
    from sklearn.compose import ColumnTransformer
    from sklearn.pipeline import make_pipeline, Pipeline
    from sklearn.ensemble import RandomForestRegressor, ExtraTreesRegressor
    from sklearn.feature_selection import RFECV, VarianceThreshold, RFE, SelectFromModel
    import scipy
    from sklearn import clone

    from sklearn.svm import SVC
    from sklearn.linear_model import LogisticRegressionCV, LogisticRegression
    from sklearn.metrics import log_loss, roc_auc_score
    from sklearn.multiclass import OneVsRestClassifier
    from sklearn.naive_bayes import MultinomialNB
    from sklearn.impute import SimpleImputer
    from sklearn.preprocessing import StandardScaler
    from sklearn.compose import ColumnTransformer
    from sklearn.pipeline import make_pipeline, Pipeline
    from sklearn.model_selection import cross_val_score, KFold, cross_val_predict, train_test_split
    from sklearn.preprocessing import LabelEncoder, OneHotEncoder, StandardScaler
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.preprocessing import MultiLabelBinarizer
    from sklearn.neural_network import MLPClassifier

    from tensorflow.python.keras.layers import Input, Dense, Activation, Flatten, Conv2D
    from tensorflow.python.keras.wrappers.scikit_learn import KerasRegressor
    from tensorflow.python.keras.models import Model, load_model
    from tensorflow import saved_model

    from tensorflow.python.keras.regularizers import l2, l1
    from Bio.SubsMat import MatrixInfo as matlist
    from sklearn.ensemble import GradientBoostingClassifier
    from sklearn.multioutput import MultiOutputClassifier
    import pickle
    # import lightgbm as lgb
    import reduce_memory_usage
    map_residue_variety_to_features()
