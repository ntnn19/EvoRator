
import click


@click.command()
@click.argument('features', type=click.Path(exists=True))
@click.argument('consurf_output_path', type=click.Path(exists=True))
@click.argument('results_dir', type=click.Path(exists=True))
@click.argument('pdb_id', type=str)
def map_residue_variety_to_features(features,consurf_output_path,results_dir,pdb_id):
    consurf_output = pd.read_csv(consurf_output_path,sep='\t',skiprows=15,skipfooter=4,header=None)
    consurf_output = consurf_output.dropna(axis=1)
    evorator_output = pd.read_csv(features)
    FINAL_OUTPUT = os.path.join(results_dir,pdb_id+"_residue_variety.csv")
    consurf_output = consurf_output[[0, 1, 2, 3, 5, 7, 9, 13, 14]]

    consurf_output.columns = ['position', 'seq', 'pdb_seq', 'normalized_score', 'color','confidence_intervals','confidence_colors', 'msa_data', 'residue_variety']

    # ALA117: A
    consurf_output = consurf_output[consurf_output['pdb_seq'] != "         -"]
    consurf_output['chain'] = consurf_output['pdb_seq'].str.split(":").str[-1]
    consurf_output['pdb_position_consurf'] = consurf_output['pdb_seq'].str.split(":").str[0].str.extract('(\d+)')
    consurf_output['pdb_position'] = consurf_output['chain'].astype(str) + consurf_output['pdb_position_consurf'].astype(str) + "_" + pdb_id
    # print(consurf_output.columns)
    # print(consurf_output['residue_variety'])
    consurf_output['residue_variety_parsed'] = consurf_output['residue_variety'].str.split(",").str.join("")
    # print(consurf_output['residue_variety_parsed'])
    # for c in 'ACDEFGHIKLMNPQRSTVWY':
    #     consurf_output[c] = 0

    mlb = MultiLabelBinarizer()
    aa_letters = 'ABCDEFGHIKLMNPQRSTVWXYZ'
    mlb.fit(aa_letters.split())
    residue_variety_parsed_binned= mlb.transform(consurf_output['residue_variety_parsed'])

    # consurf_output['ABCDEFGHIJKLMNOPQRSTUVWXYZ'.split()] = residue_variety_parsed_binned
    # print(consurf_output.columns)
    # print(mlb.classes_)

    y= pd.DataFrame(residue_variety_parsed_binned)
    print(y)
    y.columns = list(aa_letters)
    print(y)
    print(consurf_output)
    consurf_output_with_binned_labels = pd.concat([consurf_output.reset_index(),y.reset_index()],axis=1)
    print(consurf_output_with_binned_labels)
    consurf_output_with_binned_labels = consurf_output_with_binned_labels.dropna(axis=1)
    print(consurf_output_with_binned_labels)
    merged_df = evorator_output.merge(consurf_output_with_binned_labels[['pdb_position','normalized_score',*list(aa_letters)]],on='pdb_position',how='left')
    merged_df.to_csv(FINAL_OUTPUT,index=False)
    #
    # y_cols = list('ABCDEFGHIJKLMNOPQRSTUVWXYZ')
    # y = merged_df[y_cols].values
    # not_features = ['index', 'pdb_position', 'chain', 'pdb', 'pdb_x', 'pdb_y', 'chain_x', 'chain_y', 'pdb_aa_y','Functional', 'd',
    #                 '1/d','A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M',
    #                 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X','Y', 'Z','pdb_aa_x', 'total_nan_neigh', 'zr4s_JC'] + [f for f in merged_df.columns.tolist() if f.startswith('neighbor_pos_') or f.startswith("n2v")]
    #
    # categorical_features = ['glycosylation', 'protein_interaction', 'disorder', 'binding', 'catalysis', 'pdb_aa',
    #                         'aa_group_5', 'aa_group_HP', 'structure']
    #
    # all_features = [f for f in merged_df.columns.tolist() if f not in not_features]
    # numeric_features = [f for f in all_features if f not in categorical_features]
    #
    # X = merged_df[[*numeric_features,*categorical_features]]
    # print(X.columns.tolist())
    # ohe = OneHotEncoder()
    # numeric_transformer_rf = Pipeline(steps=[('imputer', SimpleImputer(strategy='median')),('scaler', StandardScaler())])
    # categorical_transformer = Pipeline(steps=[('onehot', ohe)])
    # preprocessor_rf = ColumnTransformer(
    #             transformers=[
    #                 ('num', numeric_transformer_rf, numeric_features),
    #                 ('cat', categorical_transformer, categorical_features),
    #             ])
    # steps = [('preprocessor', preprocessor_rf),('clf', OneVsRestClassifier(RandomForestClassifier(n_estimators=1000)))]
    # pipeline = Pipeline(steps)
    # X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.05, random_state=42)
    # pipeline.fit(X_train,y_train)
    # y_hat = pipeline.predict_proba(X_test)
    # print(y_hat)
    # print()
    # print(y_test)
    # S = log_loss(y_test,y_hat)
    # print(S)
    #
    # S = roc_auc_score(y_test.ravel(),y_hat.ravel())
    # print(S)
    # # S = cross_val_predict(pipeline,X,y,cv=KFold(10,shuffle=True),method='predict_proba')
    # true_scores = []
    # shuff_scores = []
    # F = True
    # if F == True:
    #     for i in range(10):
    #         try:
    #             y_shuff = y.copy()
    #             np.random.shuffle(y_shuff)
    #             X_train, X_test, y_train, y_test = train_test_split(X, y_shuff, test_size=0.05, random_state=None)
    #             pipeline.fit(X_train, y_train)
    #             y_hat = pipeline.predict_proba(X_test)
    #             S = log_loss(y_test, y_hat)
    #             S = roc_auc_score(y_test.ravel(), y_hat.ravel())
    #             shuff_scores.append(S)
    #             print(S)
    #         except:
    #             continue
    #
    #     for i in range(10):
    #         try:
    #             X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.05, random_state=None)
    #             pipeline.fit(X_train, y_train)
    #             y_hat = pipeline.predict_proba(X_test)
    #             S = log_loss(y_test, y_hat)
    #             S = roc_auc_score(y_test.ravel(), y_hat.ravel())
    #             true_scores.append(S)
    #             print(S)
    #         except:
    #             continue
    #     print(true_scores,shuff_scores)
    #
    #
    # y_hat_df = pd.DataFrame(y_hat,columns=y_cols)
    # y_test_df = pd.DataFrame(y_test,columns=y_cols)
    # d_df= pd.DataFrame(y_test_df.values-y_hat_df.values,columns=y_cols)
    # print(d_df)
    # exit()
    # d_df.to_csv("predict_aa.csv")
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
        KFold, GridSearchCV, GroupShuffleSplit
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
    from sklearn.externals import joblib

    from sklearn.svm import SVC
    from sklearn.linear_model import LogisticRegressionCV
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
    map_residue_variety_to_features()
