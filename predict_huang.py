import logging

import click


@click.command()
@click.argument('huang_features', type=click.Path(exists=True))
@click.argument('huang_rates', type=click.Path(exists=True))
@click.argument('model', type=click.Path(exists=True))
@click.argument('outdir', type=click.Path(exists=True))
def predict_huang(huang_features,huang_rates,model,outdir):

    logging.basicConfig(filename=outdir + '/' + 'huang_prediction.log',
                        level=logging.DEBUG,
                        format='%(asctime)s %(levelname)-8s %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')

    huang_features = pd.read_csv(huang_features)
    loaded_model = joblib.load(model)
    logging.debug(f'model parameters = {loaded_model.get_params()}')
    rates = pd.read_csv(huang_rates)
    rates['pdb_position'] = rates.chain.astype(str) + rates.site.astype(str) + "_" + rates.pdb.astype(str)

    FINAL_OUTPUT = os.path.join(outdir,'huang_with_evorator_features_and_predictions.csv')


    not_features = ['pdb_position', 'chain_x', 'chain', 'chain_y', 'pdb_aa_x', 'pdb_aa_y'] + [f for f in huang_features.columns.tolist() if f.startswith('neighbor_pos_')] + [f for f in huang_features.columns.tolist() if f.startswith('n2v')]
    # categorical_features = ['glycosylation', 'protein_interaction', 'disorder', 'binding', 'catalysis', 'pdb_aa_x', 'aa_group_5', 'aa_group_HP','structure']
    categorical_features = ['glycosylation', 'protein_interaction', 'disorder', 'binding', 'catalysis', 'pdb_aa',
                            'aa_group_5', 'aa_group_HP', 'structure']
    all_features = [f for f in huang_features.columns.tolist() if f not in not_features]
    numeric_features = [f for f in all_features if f not in categorical_features]
    huang_features['pdb_aa'] = huang_features['pdb_aa'].fillna(huang_features['pdb_aa_x'].fillna(huang_features['pdb_aa_y']))
    logging.debug(huang_features['pdb_aa'])
    for cat in categorical_features:
        huang_features[cat].fillna('missing', inplace=True)

    logging.debug(f'numeric\n\n{numeric_features}')
    logging.debug(f'cat\n\n{categorical_features}')
    all_features = [*numeric_features, *categorical_features]
    complete_features = ['2_clique', '3_clique', '4_clique', '5_clique', '6_clique', '7_clique',
                         'average_neighbor_degree',
                         'betweenness_centrality', 'clustering_coefficient', 'degree_centrality',
                         'eigenvector_centrality',
                         'graphlet1',
                         'graphlet10', 'graphlet11', 'graphlet12', 'graphlet13', 'graphlet14', 'graphlet15',
                         'graphlet16',
                         'graphlet17',
                         'graphlet18', 'graphlet19', 'graphlet2', 'graphlet20', 'graphlet21', 'graphlet22',
                         'graphlet23',
                         'graphlet24',
                         'graphlet25', 'graphlet26', 'graphlet27', 'graphlet28', 'graphlet29', 'graphlet3',
                         'graphlet30',
                         'graphlet31',
                         'graphlet32', 'graphlet33', 'graphlet34', 'graphlet35', 'graphlet36', 'graphlet37',
                         'graphlet38',
                         'graphlet39',
                         'graphlet4', 'graphlet40', 'graphlet41', 'graphlet42', 'graphlet43', 'graphlet44',
                         'graphlet45',
                         'graphlet46',
                         'graphlet47', 'graphlet48', 'graphlet49', 'graphlet5', 'graphlet50', 'graphlet51',
                         'graphlet52',
                         'graphlet53',
                         'graphlet54', 'graphlet55', 'graphlet56', 'graphlet57', 'graphlet58', 'graphlet59',
                         'graphlet6',
                         'graphlet60',
                         'graphlet61', 'graphlet62', 'graphlet63', 'graphlet64', 'graphlet65', 'graphlet66',
                         'graphlet67',
                         'graphlet68',
                         'graphlet69', 'graphlet7', 'graphlet70', 'graphlet71', 'graphlet72', 'graphlet73', 'graphlet8',
                         'graphlet9',
                         'median_rsa_neigh', 'median_wcn_ca_neigh', 'median_wcn_sc_neigh', 'node_degree', 'rsa',
                         'total_A_neigh',
                         'total_Aliphatic_neigh', 'total_Aromatic_neigh', 'total_C_neigh', 'total_Charged_neigh',
                         'total_D_neigh',
                         'total_Diverse_neigh', 'total_E_neigh', 'total_F_neigh', 'total_G_neigh', 'total_H_neigh',
                         'total_Hydrophobic_neigh', 'total_I_neigh', 'total_K_neigh', 'total_L_neigh', 'total_M_neigh',
                         'total_N_neigh',
                         'total_P_neigh', 'total_Polar_neigh', 'total_Q_neigh', 'total_R_neigh', 'total_S_neigh',
                         'total_T_neigh',
                         'total_Tiny_neigh', 'total_V_neigh', 'total_W_neigh', 'total_Y_neigh', 'total_catalytic_neigh',
                         'total_contact_neigh', 'total_disordered_neigh', 'total_glycosylated_neigh',
                         'total_interface_neigh',
                         'total_nan_neigh', 'total_site_neigh', 'total_ss_B_neigh', 'total_ss_E_neigh',
                         'total_ss_G_neigh',
                         'total_ss_H_neigh', 'total_ss_I_neigh', 'total_ss_P_neigh', 'total_ss_S_neigh',
                         'total_ss_T_neigh',
                         'total_ss_nan_neigh', 'wcn_ca', 'wcn_sc', 'glycosylation', 'protein_interaction', 'disorder',
                         'binding',
                         'catalysis', 'pdb_aa', 'aa_group_5', 'aa_group_HP', 'structure']


    rates_d = dict(zip(rates['pdb_position'],rates['zr4s_JC']))
    r2_scores=[]
    predictions_d = []
    i=0
    for p in huang_features.pdb.unique().tolist():
        # if i==3:
        #     break
        logging.debug(f'PDB {p}')
        curr_prot = huang_features[huang_features.pdb==p]
        X_test = curr_prot[all_features]
        curr_prot['y_hat'] = loaded_model.predict(X_test)
        curr_prot['zr4s_JC'] = curr_prot['pdb_position'].map(rates_d)
        curr_prot_clean = curr_prot[['zr4s_JC','y_hat']].dropna()
        score = r2_score(curr_prot_clean['zr4s_JC'],curr_prot_clean['y_hat'])
        logging.debug(f'r2 score = {score}')
        r2_scores.append(score)
        predictions_d.append(dict(zip(curr_prot['pdb_position'],curr_prot['y_hat'])))
        i+=1

    logging.debug(f'mean r2 = {np.mean(r2_scores)}')
    logging.debug(f'std r2 = {np.std(r2_scores)}')
    predictions_d = dict(ChainMap(*predictions_d))
    huang_features['predicted_score'] = huang_features.pdb_position.map(predictions_d)
    logging.debug(f'mean r2 = {np.mean(r2_scores)}')

    huang_features.to_csv(FINAL_OUTPUT,index=False)

if __name__ == '__main__':
    import itertools
    import time
    import logging
    import numpy as np
    import seaborn as sns
    import pandas as pd
    import pickle
    from sklearn.externals import joblib
    import os
    from scipy.stats import loguniform
    from sklearn.metrics import r2_score
    from collections import ChainMap

    predict_huang()
