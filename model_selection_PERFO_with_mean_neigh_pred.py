def r2_per_group(data, truth, predicted):
    return r2_score(data[truth], data[predicted])


def mean_confidence_interval(data, confidence=0.95):
    a = 1.0 * np.array(data)
    n = len(a)
    m, se = np.mean(a), scipy.stats.sem(a)
    h = se * scipy.stats.t.ppf((1 + confidence) / 2., n - 1)
    return m, m - h, m + h


import click


@click.command()
@click.argument('features', type=click.Path(exists=True))
@click.argument('huang_rates', type=click.Path(exists=True))
@click.argument('results_dir', type=click.Path(exists=True))
@click.option('--use-traditional-features', is_flag=True, default=True,
              help='If given uses all features except 4 graph features')
@click.option('--pred', is_flag=True, default=False, help='If given predicts with selected model')
@click.option('--fit', is_flag=True, default=False, help='If given fits model to data')
@click.option('--trained-model', type=click.Path(exists=True), help='path to a trained model')
@click.option('--estimator',
              type=click.Choice(['RF', 'SVR', 'LinearRegression', 'Lasso'], case_sensitive=False))
@click.option('--lasso-regularization', type=float)
@click.option('--rf-n-trees', type=int)
@click.option('--rf-max-depth', type=int)
@click.option('--svr-c', type=float)
@click.option('--svr-gamma', type=click.Choice(['scale', '0.00001', '0.0001', '0.001', '0.01', '0.1', '1', '10']))
def main(features, huang_rates, use_traditional_features, results_dir, estimator, lasso_regularization, rf_n_trees,
         rf_max_depth, svr_gamma, svr_c, fit, pred, trained_model):
    ''' Uses the hold-out set method to evaluate the preformance of traditional vs evorator feature set'''

    # prepare dataset
    if svr_gamma:
        try:
            svr_gamma = float(svr_gamma)
        except:
            pass

    if use_traditional_features:
        suffix = "all_features"
    else:
        suffix = "traditional_features"
    features = pd.read_csv(features)
    rates = pd.read_csv(huang_rates)
    # selected_protein = prot_d.get(test_protein)
    # selected_protein_dir = os.path.join(results_dir,selected_protein)
    # os.makedirs(selected_protein_dir,exist_ok=True)
    logging.basicConfig(
        filename=results_dir + '/' + suffix + "_neigh_pred_" + estimator + '_model_selection.log',
        level=logging.DEBUG,
        format='%(asctime)s %(levelname)-8s %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S')

    # logging.debug(f'{prot_d}')
    # logging.debug(f'selected: {selected_protein}')
    # df_merged = df_merged[df_merged['chain_x']==pdb_chain]

    # A2_132L

    rates['pdb_position'] = rates.chain.astype(str) + rates.site.astype(str) + "_" + rates.pdb.astype(str)
    training_set = rates.drop(columns=['chain', 'site', 'pdb', 'zr4s_JTT', 'r4s_JC', 'r4s_JTT']).merge(features,
                                                                                                       on='pdb_position')
    logging.debug(training_set)
    logging.debug(training_set.pdb.unique())
    logging.debug(training_set.pdb.unique().shape)
    training_set['chain'] = training_set['chain'].fillna(training_set['chain_x'])
    training_set['chain'] = training_set['chain'].fillna(training_set['chain_y'])
    training_set['pdb_aa'] = training_set['pdb_aa'].fillna(training_set['pdb_aa_x'])
    training_set['pdb_aa'] = training_set['pdb_aa'].fillna(training_set['pdb_aa_y'])

    # logging.debug(f'MERGED TABLE COLUMNS: {training_set.columns.tolist()}')
    logging.debug(f'MERGED TABLE: {training_set.columns.tolist()}')
    not_features = ['index', 'pdb_position', 'chain', 'pdb', 'pdb_x', 'pdb_y', 'chain_x', 'chain_y', 'pdb_aa_y',
                    'pdb_aa_x', 'total_nan_neigh', 'zr4s_JC'] + [f for f
                                                                 in
                                                                 training_set.columns.tolist()
                                                                 if
                                                                 f.startswith(
                                                                     'neighbor_pos_') or f.startswith("n2v")] \
        # + [f for f
    #                                               in
    #                                               training_set.columns.tolist()
    #                                               if
    #                                               f.endswith(
    #                                                   '_neigh') and f.startswith("total")]

    categorical_features = ['glycosylation', 'protein_interaction', 'disorder', 'binding', 'catalysis', 'pdb_aa',
                            'aa_group_5', 'aa_group_HP', 'structure']

    all_features = [f for f in training_set.columns.tolist() if f not in not_features]
    numeric_features = [f for f in all_features if f not in categorical_features]


    if not use_traditional_features:
        graph_features = [f for f in numeric_features if "graphlet" in f]
        graph_features = [*graph_features, *[f for f in numeric_features if "clique" in f]]
        graph_features = [*graph_features, *[f for f in numeric_features if "_neigh" in f]]
        graph_features = [*graph_features, 'node_degree', 'average_neighbor_degree', 'betweenness_centrality',
                          'clustering_coefficient', 'degree_centrality', 'eigenvector_centrality']

        n2v_features = [f for f in numeric_features if "n2v" in f]

        numeric_features = [f for f in all_features if f not in [*categorical_features, *graph_features, *n2v_features]]
        categorical_features = []  # remove this line to include more traditional features
        numeric_features = ['rsa', 'wcn_ca']  # remove this line to include more numerical features

    all_features = [*numeric_features, *categorical_features]
    # logging.debug(f'selected features:\n\n {all_features}')
    logging.debug(f'features:\n\n {all_features}')
    if use_traditional_features:  # remove condition to include more traditional features
        for cat in categorical_features:
            training_set[cat].fillna('missing', inplace=True)
    logging.debug(f'cat features:{categorical_features}')
    logging.debug(f'numeric_features : {numeric_features}')
    training_set_clean = training_set.dropna(subset=['zr4s_JC']).reset_index()
    ohe = OneHotEncoder(handle_unknown='ignore')
    numeric_transformer = Pipeline(steps=[('imputer', SimpleImputer(strategy='median')), ('scaler', StandardScaler())])
    numeric_transformer_rf = Pipeline(steps=[('imputer', SimpleImputer(strategy='median'))])

    if use_traditional_features:  # remove condition to include more traditional features
        categorical_transformer = Pipeline(steps=[('onehot', ohe)])

    if use_traditional_features:  # remove condition to include more traditional features
        preprocessor = ColumnTransformer(
            transformers=[
                ('num', numeric_transformer, numeric_features),
                ('cat', categorical_transformer, categorical_features),
            ])
        preprocessor_rf = ColumnTransformer(
            transformers=[
                ('num', numeric_transformer_rf, numeric_features),
                ('cat', categorical_transformer, categorical_features),
            ])
    else:
        preprocessor = ColumnTransformer(
            transformers=[
                ('num', numeric_transformer, numeric_features),
            ])
        preprocessor_rf = ColumnTransformer(
            transformers=[
                ('num', numeric_transformer_rf, numeric_features),
            ])
    if estimator == 'RF':
        regression_algorithm = RandomForestRegressor(n_estimators=rf_n_trees, max_depth=rf_max_depth)
        preprocessor = preprocessor_rf
    elif estimator == 'SVR':
        regression_algorithm = SVR(C=svr_c, gamma=svr_gamma)
    elif estimator == 'LinearRegression':
        regression_algorithm = LinearRegression()
    elif estimator == 'Lasso':
        regression_algorithm = Lasso(alpha=lasso_regularization, max_iter=5000)

    X = training_set_clean[[*all_features]]
    y = training_set_clean['zr4s_JC']
    logging.debug(f'{X.shape} {y.shape}')

    prot_d = dict(zip(training_set_clean['pdb'].unique().tolist(), [str(i) for i in range(1, 214)]))

    import random
    random.seed(2389)
    rand_idx = list(random.sample(range(213), 163))
    train_proteins = training_set_clean['pdb'].unique()[rand_idx]
    rand_idx = list(random.sample(range(163), 30))
    val_proteins = train_proteins[rand_idx]
    train_proteins = [p for p in train_proteins if p not in val_proteins]
    test_proteins = [p for p in training_set_clean['pdb'].unique() if p not in [*train_proteins, *val_proteins]]

    train_ix = training_set_clean.index[training_set_clean['pdb'].isin(train_proteins)].tolist()
    test_ix = training_set_clean.index[training_set_clean['pdb'].isin(test_proteins)].tolist()
    val_ix = training_set_clean.index[training_set_clean['pdb'].isin(val_proteins)].tolist()
    logging.debug(
        f'train proteins={training_set_clean["pdb"].iloc[train_ix].unique()},len={training_set_clean["pdb"].iloc[train_ix].unique().shape}')  # 0.25 x 0.8 = 0.2
    logging.debug(
        f'val proteins={training_set_clean["pdb"].iloc[val_ix].unique()},len={training_set_clean["pdb"].iloc[val_ix].unique().shape}')  # 0.25 x 0.8 = 0.2
    logging.debug(
        f'test proteins={training_set_clean["pdb"].iloc[test_ix].unique()},len={training_set_clean["pdb"].iloc[test_ix].unique().shape}')  # 0.25 x 0.8 = 0.2

    X_train = X.iloc[train_ix, :]
    y_train = y.iloc[train_ix]
    X_test = X.iloc[test_ix, :]
    y_test = y.iloc[test_ix]
    X_val = X.iloc[val_ix, :]
    y_val = y.iloc[val_ix]

    # logging.debug(sum([1 if v in train_ix else 0 for v in test_ix]))
    # logging.debug(sum([1 if v in train_ix else 0 for v in val_ix]))
    logging.debug(X_train.shape)  # 0.25 x 0.8 = 0.2
    logging.debug(X_val.shape)  # 0.25 x 0.8 = 0.2
    logging.debug(X_test.shape)  # 0.25 x 0.8 = 0.2
    # Data Prep
    preprocessor.fit(X_train)
    X_train_enc = preprocessor.transform(X_train)
    X_val_enc = preprocessor.transform(X_val)

    # feature selection
    t1 = time.perf_counter()
    if estimator == 'SVR':
        fs_algorithm = SelectFromModel(LassoCV(
            alphas=[10 ** -8, 10 ** -7, 10 ** -6, 10 ** -5, 10 ** -4, 10 ** -3, 10 ** -2, 10 ** -1, 10 ** 0, 10 ** 1,
                    10 ** 2], random_state=216, max_iter=10 ** 5))

        if fit and pred:
            steps = [('preprocessor', preprocessor),
                     ('FS', fs_algorithm),
                     ('SVR', SVR(C=svr_c, gamma=svr_gamma))]
            pipeline = Pipeline(steps)


            # X_train = pd.concat([X_train,X_val]).iloc[:1000,:]
            X_train = pd.concat([X_train,X_val])
            # y_train = pd.concat([y_train,y_val]).iloc[:1000]
            y_train = pd.concat([y_train,y_val])
            X_test = X.iloc[test_ix, :]
            y_test = y.iloc[test_ix]

            logging.debug(f'{X_train.columns.tolist()}')
            logging.debug(f'fitting 231')
            pipeline.fit(X_train, y_train)
            logging.debug(f'predicting')
            y_hat = pipeline.predict(X_test)
            training_score = r2_score(y_test, y_hat)
            logging.debug(f'test score={training_score}')

            # training_set = pd.concat([training_set_clean.iloc[train_ix, :].copy() ,training_set_clean.iloc[val_ix, :].copy()],axis=0).iloc[:1000,:]
            training_set = pd.concat([training_set_clean.iloc[train_ix, :].copy() ,training_set_clean.iloc[val_ix, :].copy()],axis=0)
            test_set = training_set_clean.iloc[test_ix, :].copy()



            training_set['y_hat'] = pipeline.predict(X_train)
            test_set['y_hat'] = y_hat

            # Add feature to training data
            training_neigh_df = training_set[['pdb_position', *[f for f in training_set.columns.tolist() if 'neighbor_pos_' in f], 'y_hat']]
            training_y_hat_d = dict(zip(training_neigh_df['pdb_position'], training_neigh_df['y_hat']))
            training_neigh_d = dict(zip(training_neigh_df['pdb_position'],training_neigh_df[[f for f in training_set.columns.tolist() if 'neighbor_pos_' in f]].values.tolist()))
            training_mean_neigh_y_hat_d = {}
            for k in training_neigh_d:
                tmp = []
                for n in training_neigh_d[k]:
                    tmp.append(training_y_hat_d.get(n, np.nan))
                training_mean_neigh_y_hat_d[k] = np.nanmean(tmp)
            training_set['mean_neigh_y_hat'] = training_set['pdb_position'].map(training_mean_neigh_y_hat_d)

            # Add feature to test data
            test_neigh_df = test_set[['pdb_position', *[f for f in test_set.columns.tolist() if 'neighbor_pos_' in f], 'y_hat']]
            test_y_hat_d = dict(zip(test_neigh_df['pdb_position'], test_neigh_df['y_hat']))
            test_neigh_d = dict(zip(test_neigh_df['pdb_position'], test_neigh_df[
                [f for f in test_set.columns.tolist() if 'neighbor_pos_' in f]].values.tolist()))
            test_mean_neigh_y_hat_d = {}
            for k in test_neigh_d:
                tmp = []
                for n in test_neigh_d[k]:
                    tmp.append(test_y_hat_d.get(n, np.nan))
                test_mean_neigh_y_hat_d[k] = np.nanmean(tmp)
            test_set['mean_neigh_y_hat'] = test_set['pdb_position'].map(test_mean_neigh_y_hat_d)

            # Evaluate mean predicted score contribution
            preprocessor_mean_neigh_yhat = ColumnTransformer(
                transformers=[
                    ('num', numeric_transformer, ['mean_neigh_y_hat','y_hat']),
                ])
            steps_mean_neigh_yhat = [('preprocessor', preprocessor_mean_neigh_yhat),
                     ('SVR', SVR())]


            pipeline_mean_neigh_yhat = Pipeline(steps_mean_neigh_yhat)

            X_train_mean_neigh_yhat = training_set[['mean_neigh_y_hat','y_hat']]
            y_train_mean_neigh_yhat = training_set['zr4s_JC']
            # X_val_mean_neigh_yhat = training_set[['mean_neigh_y_hat','y_hat']].iloc[val_ix,:]
            # y_val_mean_neigh_yhat = training_set['zr4s_JC'].iloc[val_ix]
            # X_train_mean_neigh_yhat = pd.concat([X_train_mean_neigh_yhat, X_val_mean_neigh_yhat]).iloc[:1000]
            # y_train_mean_neigh_yhat = pd.concat([y_train_mean_neigh_yhat, y_val_mean_neigh_yhat]).iloc[:1000]
            X_test_mean_neigh_yhat = test_set[['mean_neigh_y_hat','y_hat']]
            y_test_mean_neigh_yhat = test_set['zr4s_JC']

            pipeline_mean_neigh_yhat.fit(X_train_mean_neigh_yhat,y_train_mean_neigh_yhat)
            y_hat_test_mean_neigh_yhat = pipeline_mean_neigh_yhat.predict(X_test_mean_neigh_yhat)
            test_score_mean_neigh_yhat = r2_score(y_test_mean_neigh_yhat, y_hat_test_mean_neigh_yhat)
            y_hat_train_mean_neigh_yhat = pipeline_mean_neigh_yhat.predict(X_train_mean_neigh_yhat)
            training_score_mean_neigh_yhat = r2_score(y_train_mean_neigh_yhat, y_hat_train_mean_neigh_yhat)

            logging.debug(f'training_mean_neigh_yhat score ={training_score_mean_neigh_yhat}')
            logging.debug(f'test_mean_neigh_yhat score ={test_score_mean_neigh_yhat}')

            logging.debug(f'saving features & predictions')
            test_set['new_pred'] = y_hat_test_mean_neigh_yhat
            training_set['new_pred'] = y_hat_train_mean_neigh_yhat
            test_set.to_csv(os.path.join(results_dir, f'EvoRator_SVR_C=default_gamma=default_model_PERFO_test_set_pred_mean_neigh_yhat_final.csv'), index=False)
            training_set.to_csv(os.path.join(results_dir, f'EvoRator_SVR_C=default_gamma=default_model_PERFO_training_set_pred_mean_neigh_yhat_final.csv'), index=False)
            F = False
            if F:
                logging.debug(f'saving model')
                joblib.dump(pipeline, os.path.join(results_dir, f'EvoRator_SVR-FS_C={svr_c},gamma={svr_gamma}_model_PERFO' + '.joblib'))
            logging.debug(f'done')
            return
        elif pred and trained_model:

            logging.debug(f'Just predicting')
            logging.debug(f'loading model')
            pipeline = joblib.load(trained_model)
            logging.debug(f'predicting')
            y_hat = pipeline.predict(X)
            training_score = r2_score(y, y_hat)
            logging.debug(f'training score={training_score}')
            training_set_clean['y_hat'] = y_hat
            logging.debug(f'saving features & predictions')
            training_set_clean.to_csv(os.path.join(results_dir, f'EvoRator_SVR-FS_C={svr_c},gamma={svr_gamma}_model_PERFGR_pred' + '.csv'))
            return



    # evaluate the model
    # fit the model
    logging.debug(f'Feature selection')
    fs_algorithm.fit(X_train_enc, y_train)
    X_train_enc = fs_algorithm.transform(X_train_enc)

    skip_block = True
    if not skip_block:
        X_val_enc = fs_algorithm.transform(X_val_enc)
        logging.debug(f'fitting')
        regression_algorithm.fit(X_train_enc, y_train)
        # train score
        yhat_train = regression_algorithm.predict(X_train_enc)
        training_score = r2_score(y_train, yhat_train)
        # validation score
        logging.debug(f'predicting using the validation set')
        yhat_val = regression_algorithm.predict(X_val_enc)
        val_score = r2_score(y_val, yhat_val)
        logging.debug(f'training_score= {training_score} , validation_score= {val_score}')

    if use_traditional_features:
        categorical_feature_names = preprocessor.transformers_[1][1].named_steps['onehot'].get_feature_names(
            categorical_features)
        selected_features = list(np.array(list(numeric_features) + list(categorical_feature_names))[fs_algorithm.get_support()])
        logging.debug(
            f"The selected features are {list(np.array(list(numeric_features) + list(categorical_feature_names))[fs_algorithm.get_support()])}")
    else:

        logging.debug(
            f"The selected features are {list(np.array(list(numeric_features))[fs_algorithm.get_support()])}")
# return
# Use the selected parameters to retrain the algorithm on all training data and evaluate on the test set

    X_train = pd.concat([X_train, X_val], axis=0).iloc[:1000,:]
    y_train = pd.concat([y_train, y_val], axis=0).iloc[:1000]
    # X_train = pd.concat([X_train, X_val], axis=0)
    # y_train = pd.concat([y_train, y_val], axis=0)
    preprocessor.fit(X_train)
    X_train_enc = preprocessor.transform(X_train)
    X_test_enc = preprocessor.transform(X_test)
    logging.debug(f'Use the selected features')
    X_train_enc = X_train_enc[selected_features]
    X_test_enc = X_test_enc[selected_features]


    # train score
    regression_algorithm = clone(regression_algorithm)
    logging.debug(f'fitting')
    regression_algorithm.fit(X_train_enc, y_train)
    yhat_train = regression_algorithm.predict(X_train_enc)
    training_score = r2_score(y_train, yhat_train)

    # test score
    logging.debug(f'predicting')
    yhat_test = regression_algorithm.predict(X_test_enc)
    test_score = r2_score(y_test, yhat_test)

    # report progress
    t2 = time.perf_counter()
    logging.debug(f'training_score= {training_score}')
    logging.debug(f'test_score= {test_score}')
    logging.debug(f'Run time={t2 - t1}')

    df_test = training_set_clean.iloc[test_ix, :].copy()
    df_test['pred'] = yhat_test
    logging.debug(f'saving features and predictions of the hold-out set')
    results_path = os.path.join(results_dir, suffix + f"_test_set_pred_{estimator}_for2step_final.csv")
    df_test.to_csv(results_path, index=False)
    return # remove this line to write model selection results
    dfs = []
    for name, t in zip(['test'], [df_test]):
        tg = t[['pdb', 'zr4s_JC', 'pred']].groupby('pdb')
        n = tg.apply(r2_per_group, 'zr4s_JC', 'pred')
        n.name = name
        logging.debug(f'{n}')
        # res = bootstrap((n,), np.median)
        logging.debug(f'{name}_score= {np.median(n)}')
        logging.debug(f'{name}_iqr= {n.describe()}')

    if estimator == 'RF':
        results_df = pd.DataFrame(data=[[rf_n_trees, rf_max_depth, training_score, val_score, test_score]],
                                  columns=['n_trees', 'max_depth', 'training_score', 'val_score', 'test_score'])
        results_path = os.path.join(results_dir, suffix + "_" + str(rf_n_trees) + "_" + str(
            rf_max_depth) + f"_scores_{estimator}.csv")
        results_df_test_proteins_path = os.path.join(results_dir, suffix + "_" + str(rf_n_trees) + "_" + str(
            rf_max_depth) + f"_scores_for_proteins_{estimator}.csv")
        results_df.to_csv(results_path, index=False)
        n.to_csv(results_df_test_proteins_path)
    elif estimator == 'SVR':
        results_df = pd.DataFrame(data=[[svr_c, svr_gamma, training_score, val_score, test_score]],
                                  columns=['svr_c', 'svr_gamma', 'training_score', 'val_score', 'test_score'])
        results_path = os.path.join(results_dir,
                                    suffix + "_" + str(svr_c) + "_" + str(svr_gamma) + f"_scores_{estimator}.csv")
        results_df_test_proteins_path = os.path.join(results_dir, suffix + "_" + str(svr_c) + "_" + str(
            svr_gamma) + f"_scores_for_proteins_{estimator}.csv")
        results_df.to_csv(results_path, index=False)
        n.to_csv(results_df_test_proteins_path)
    elif estimator == 'Lasso':
        results_df = pd.DataFrame(data=[[lasso_regularization, training_score, val_score, test_score]],
                                  columns=['lasso_regularization', 'training_score', 'val_score', 'test_score'])
        results_df_test_proteins_path = os.path.join(results_dir, suffix + "_" + str(
            lasso_regularization) + f"_scores_for_proteins_{estimator}.csv")
        results_path = os.path.join(results_dir, suffix + "_" + str(lasso_regularization) + f"_scores_{estimator}.csv")
        results_df.to_csv(results_path, index=False)
        n.to_csv(results_df_test_proteins_path)
    elif estimator == 'LinearRegression':
        results_df = pd.DataFrame(data=[[training_score, val_score, test_score]],
                                  columns=['training_score', 'val_score', 'test_score'])
        results_df_test_proteins_path = os.path.join(results_dir, suffix + f"_scores_for_proteins_{estimator}.csv")
        results_path = os.path.join(results_dir, suffix + f"_scores_{estimator}.csv")
        results_df.to_csv(results_path, index=False)
        n.to_csv(results_df_test_proteins_path)


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

    main()
