from pytorch_tabnet.metrics import Metric
def rmspe(y_true, y_pred):
    # Function to calculate the root mean squared percentage error
    return np.sqrt(np.mean(np.square((y_true - y_pred) / y_true)))



class MAE(Metric):
    def __init__(self):
        self._name = "mae"
        self._maximize = False

    def __call__(self, y_true, y_score):
        output_errors = np.average(np.abs(y_score - y_true))
        return output_errors


def MSE(y_pred, y_true):
    output_errors = torch.mean((y_true - y_pred) ** 2)
    return output_errors.clone()

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
@click.option('--estimator', type=click.Choice(['TabNet'], case_sensitive=False))
def main(features, huang_rates, use_traditional_features, results_dir, estimator):
    ''' Uses the hold-out set method to evaluate the preformance of traditional vs evorator feature set'''

    # prepare dataset

    if use_traditional_features:
        suffix = "all_features"
    else:
        suffix = "traditional_features"

    print('loading data')
    features =  dt.fread(features).to_pandas()
    rates = dt.fread(huang_rates).to_pandas()
    print('data loaded')
    # selected_protein = prot_d.get(test_protein)
    # selected_protein_dir = os.path.join(results_dir,selected_protein)
    # os.makedirs(selected_protein_dir,exist_ok=True)
    logging.basicConfig(
        filename=results_dir + '/' + suffix + "_" + estimator + '_model_selection.log',
        level=logging.DEBUG,
        format='%(asctime)s %(levelname)-8s %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S')

    # logging.debug(f'{prot_d}')
    # logging.debug(f'selected: {selected_protein}')
    # df_merged = df_merged[df_merged['chain_x']==pdb_chain]

    # A2_132L
    flg = True
    if flg:
        f = r"C:\Users\natan\Documents\deg_project\consurf_n2v_huang\data\wcn_profiles_209_monomers.csv"

        features = dt.fread(f).to_pandas()
        features['pdb_position'] = features.chain.astype(str) + features.site.astype(str) + "_" + features.pdb.astype(str)


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
    logging.debug(f'cat features:{categorical_features}')
    logging.debug(f'numeric_features : {numeric_features}')
    training_set_clean = training_set.dropna(subset=['zr4s_JC']).reset_index()


    X = training_set_clean[[*all_features]]
    y = training_set_clean['zr4s_JC']
    logging.debug(f'{X.shape} {y.shape}')


    import random
    random.seed(2389)
    rand_idx = list(random.sample(range(213), 163))
    train_proteins = training_set_clean['pdb'].unique()[rand_idx]
    rand_idx = list(random.sample(range(163), 30))
    val_proteins = train_proteins[rand_idx][:2]
    train_proteins = [p for p in train_proteins if p not in val_proteins][:75]
    test_proteins = [p for p in training_set_clean['pdb'].unique() if p not in [*train_proteins, *val_proteins]][:2]

    train_ix = training_set_clean.index[training_set_clean['pdb'].isin(train_proteins)].tolist()
    test_ix = training_set_clean.index[training_set_clean['pdb'].isin(test_proteins)].tolist()
    val_ix = training_set_clean.index[training_set_clean['pdb'].isin(val_proteins)].tolist()
    logging.debug(
        f'train proteins={training_set_clean["pdb"].iloc[train_ix].unique()},len={training_set_clean["pdb"].iloc[train_ix].unique().shape}')  # 0.25 x 0.8 = 0.2
    logging.debug(
        f'val proteins={training_set_clean["pdb"].iloc[val_ix].unique()},len={training_set_clean["pdb"].iloc[val_ix].unique().shape}')  # 0.25 x 0.8 = 0.2
    logging.debug(
        f'test proteins={training_set_clean["pdb"].iloc[test_ix].unique()},len={training_set_clean["pdb"].iloc[test_ix].unique().shape}')  # 0.25 x 0.8 = 0.2

    # logging.debug(sum([1 if v in train_ix else 0 for v in test_ix]))
    # logging.debug(sum([1 if v in train_ix else 0 for v in val_ix]))
    if estimator=='TabNet':
        categorical_dims = {}
        le = LabelEncoder()

        for col in categorical_features:
            print(col, X[col].nunique())
            l_enc = le
            X[col] = X[col].fillna("VV_likely")
            X[col] = l_enc.fit_transform(X[col].values)
            categorical_dims[col] = len(l_enc.classes_)

        X_train = X.iloc[train_ix, :]
        y_train = y.iloc[train_ix]
        X_test = X.iloc[test_ix, :]
        y_test = y.iloc[test_ix]
        X_val = X.iloc[val_ix, :]
        y_val = y.iloc[val_ix]

        logging.debug(X_train.shape)  # 0.25 x 0.8 = 0.2
        logging.debug(X_val.shape)  # 0.25 x 0.8 = 0.2
        logging.debug(X_test.shape)  # 0.25 x 0.8 = 0.2

        for col in numeric_features:

            X_train[col] = X_train[col].fillna(X_train[col].mean())
            X_val[col] = X_val[col].fillna(X_train[col].mean())
            X_test[col] = X_test[col].fillna(X_train[col].mean())
            scaler = StandardScaler()
            X_train[col] = scaler.fit_transform(X_train[col].values.reshape(-1, 1))
            X_val[col] = scaler.transform(X_val[col].values.reshape(-1, 1))
            X_test[col] = scaler.transform(X_test[col].values.reshape(-1, 1))
            X_train_enc = X_train
            X_val_enc = X_val
            X_test_enc = X_test


    # train score
    cat_idxs = [i for i, f in enumerate(X.columns.tolist()) if f in categorical_features]

    cat_dims = [categorical_dims[f] for i, f in enumerate(X.columns.tolist()) if f in categorical_features]

    tabnet_params = dict(
        cat_idxs=cat_idxs,
        cat_dims=cat_dims,
        cat_emb_dim=1,
        n_d=16,
        n_a=16,
        n_steps=2,
        gamma=2,
        n_independent=2,
        n_shared=2,
        lambda_sparse=0,
        optimizer_fn=Adam,
        optimizer_params=dict(lr=(2e-2)),
        mask_type="entmax",
        scheduler_params=dict(T_0=200, T_mult=1, eta_min=1e-4, last_epoch=-1, verbose=False),
        scheduler_fn=CosineAnnealingWarmRestarts,
        seed=42,
        verbose=10

    )
    logging.debug(f'fitting')

    regression_algorithm =  TabNetRegressor(**tabnet_params)
    # regression_algorithm.fit(
    #     X_train_enc.values, y_train.values.reshape(-1, 1),
    #     eval_set=[(X_val_enc.values, y_val.values.reshape(-1, 1))],
    #     max_epochs=50,
    #     patience=8,
    #     batch_size=1024 * 20,
    #     virtual_batch_size=128 * 20,
    #     num_workers=4,
    #     drop_last=False,
    #     eval_metric=[MAE],
    #     loss_fn=MSE
    # )
    regression_algorithm =  LinearRegression()
    #
    regression_algorithm.fit(X_train_enc, y_train)
    yhat_train = regression_algorithm.predict(X_train_enc.values).reshape(-1, 1)
    training_score = r2_score(y_train, yhat_train)
    print('training_score',training_score)
    training_score = rmspe(y_train, yhat_train)
    print('training_score',training_score)
    training_score = mean_squared_error(y_train, yhat_train)
    print('training_score',training_score)
    training_score = mean_absolute_error(y_train, yhat_train)
    print('training_score',training_score)

    # val score
    logging.debug(f'predicting')
    yhat_val = regression_algorithm.predict(X_val_enc).reshape(-1, 1)
    val_score = r2_score(y_val, yhat_val)
    print('val_score',val_score)
    val_score = rmspe(y_val, yhat_val)
    print('val_score',val_score)
    val_score = mean_squared_error(y_val, yhat_val)
    print('val_score',val_score)
    val_score = mean_absolute_error(y_val, yhat_val)
    print('val_score',val_score)

    # test score
    logging.debug(f'predicting')
    yhat_test = regression_algorithm.predict(X_test_enc).reshape(-1, 1)
    test_score = r2_score(y_test, yhat_test)
    print('test_score',test_score)
    test_score = rmspe(y_test, yhat_test)
    print('test_score',test_score)
    test_score = mean_squared_error(y_test, yhat_test)
    print('test_score',test_score)
    test_score = mean_absolute_error(y_test, yhat_test)
    print('test_score',test_score)
    exit()
    # report progress

    logging.debug(f'training_score= {training_score}')
    logging.debug(f'test_score= {test_score}')

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
    elif estimator == 'TabNet':
        results_df = pd.DataFrame(data=[[training_score, val_score, test_score]],
                                  columns=['training_score', 'val_score', 'test_score'])
        results_df_test_proteins_path = os.path.join(results_dir, suffix + f"_scores_for_proteins_{estimator}.csv")
        results_path = os.path.join(results_dir, suffix + f"_scores_{estimator}.csv")
        results_df.to_csv(results_path, index=False)
        n.to_csv(results_df_test_proteins_path)


if __name__ == '__main__':
    import numpy as np
    from sklearn.linear_model import LinearRegression, Lasso, Ridge, LassoCV
    import pandas as pd
    import random
    import os
    from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error
    import logging
    import pickle
    from sklearn.preprocessing import LabelEncoder, OneHotEncoder, StandardScaler
    import scipy
    import datatable as dt
    from torch.optim import Adam, SGD
    from torch.optim.lr_scheduler import ReduceLROnPlateau, CosineAnnealingWarmRestarts
    from pytorch_tabnet.tab_model import TabNetRegressor

    import torch
    main()
