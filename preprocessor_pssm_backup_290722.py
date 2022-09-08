import sys
import pickle
from sklearn.feature_selection import VarianceThreshold
from sklearn.preprocessing import LabelEncoder, OneHotEncoder, StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import make_pipeline, Pipeline
from Bio.SubsMat import MatrixInfo as matlist
import pandas as pd
import logging
import os
from tensorflow.python.keras.layers import Input, Dense, Activation, Flatten, Conv2D
from tensorflow.python.keras.wrappers.scikit_learn import KerasRegressor, KerasClassifier
from tensorflow.python.keras.models import Model, load_model

import pickle

from tensorflow.python.keras.layers import deserialize, serialize
from tensorflow.python.keras.saving import saving_utils


def unpack(model, training_config, weights):
    restored_model = deserialize(model)
    if training_config is not None:
        restored_model.compile(
            **saving_utils.compile_args_from_training_config(
                training_config
            )
        )
    restored_model.set_weights(weights)
    return restored_model

# Hotfix function
def make_keras_picklable():

    def __reduce__(self):
        model_metadata = saving_utils.model_metadata(self)
        training_config = model_metadata.get("training_config", None)
        model = serialize(self)
        weights = self.get_weights()
        return (unpack, (model, training_config, weights))

    cls = Model
    cls.__reduce__ = __reduce__


def NNModel_class(input_dim=None, num_classes=2,n_hidden=147,vote=False):
    X_input = Input(input_dim)
    X = Dense(input_dim[0],kernel_initializer='normal', activation='relu', name='fc1')(X_input)
    X = Dense(n_hidden, activation='relu', name='fc2')(X)
    # X = Dense(60, activation='relu', name='fc3')(X)
    # X = Dense(30, activation='relu', name='fc4')(X)
    # X = Dense(num_classes, activation='sigmoid', name='fc_out')(X)
    X = Dense(num_classes, activation='softmax', name='fc_out')(X) # as in preprint
    # X = Dense(num_classes,kernel_initializer='normal', name='fc_out')(X)
    if vote:
        X =Flatten()(X)

    # X = Dense(num_classes, activation='relu', name='fc_out')(X_input)
    model = Model(inputs=X_input, outputs=X)
    # model.compile('adam', 'binary_crossentropy', metrics=['binary_crossentropy'])
    model.compile('adam', 'kullback_leibler_divergence', metrics=['kullback_leibler_divergence']) # as in preprint
    # model.compile('adam', 'mean_squared_error', metrics=['mean_squared_error'])
    return model

HUANG_DATA_PATH = sys.argv[1]
HUANG_PSSM_PATH = sys.argv[2]
OUT_PATH = sys.argv[3]
OUT_NAME = sys.argv[4]
TEST_DATA_PATH = sys.argv[5]
RESULTS_DIR = sys.argv[6]
JOB_TITLE = sys.argv[7]
N_HIDDEN = 147
make_keras_picklable()

FITTED_PREPROCESSOR_OUTPATH = os.path.join(OUT_PATH,'preprocessor_pssm_'+OUT_NAME+'.pkl')
TRAINED_MODEL_OUTPATH = os.path.join(OUT_PATH,'trained_ANN_pssm_'+OUT_NAME+'.pkl')
TRAININIG_PREDICTIONS_OUTPATH = os.path.join(OUT_PATH,'huang_pssm_pred_'+OUT_NAME+'.csv')
# logging.basicConfig(filename=OUT_PATH + '/' + 'preprocessor.log',
#                     level=logging.DEBUG,
#                     format='%(asctime)s %(levelname)-8s %(message)s',
#                     datefmt='%Y-%m-%d %H:%M:%S')

if not os.path.exists(FITTED_PREPROCESSOR_OUTPATH) and not os.path.exists(TRAINED_MODEL_OUTPATH):
    test_data = pd.read_csv(TEST_DATA_PATH)

    huang_df = pd.read_csv(HUANG_DATA_PATH)
    huang_pssm_df = pd.read_csv(HUANG_PSSM_PATH)
    scores_d = dict(zip(huang_df['pdb_position'].tolist(),huang_df['predicted_score'].tolist()))
    huang_pssm_df['predicted_evorator_score'] = huang_pssm_df['pdb_position'].map(scores_d)
    matrix = matlist.blosum62
    aa = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

    blosum_d = {}
    for AA1 in aa:
        tmp_vec = []
        for AA2 in aa:
            tmp_vec.append(matrix.get((AA1, AA2), matrix.get((AA2, AA1))))
        blosum_d[AA1] = tmp_vec
    blosum_cols = ['blosum' + str(i) for i in range(len(blosum_d['A']))]
    huang_pssm_df[blosum_cols] = huang_pssm_df['pdb_aa'].map(blosum_d).apply(pd.Series)

    aa_freqs = {'A': 8.25 / 100, 'Q': 3.93 / 100, 'L': 9.65 / 100, 'S': 6.63 / 100, 'R': 5.53 / 100,
                'E': 6.72 / 100, 'K': 5.80 / 100, 'T': 5.35 / 100, 'N': 4.05 / 100,
                'G': 7.08 / 100, 'M': 2.41 / 100, 'W': 1.09 / 100, 'D': 5.46 / 100, 'H': 2.27 / 100,
                'F': 3.86 / 100, 'Y': 2.92 / 100, 'C': 1.38 / 100, 'I': 5.91 / 100, 'P': 4.73 / 100,
                'V': 6.86 / 100}
    huang_pssm_df['aa_freq_ExPASy'] = huang_pssm_df['pdb_aa'].map(aa_freqs)
    aa = [c + "_freq" for c in aa]
    not_features = ['index', 'pdb_position', 'pdb_id', 'chain', 'pdb', 'pdb_x', 'pdb_y', 'chain_x', 'chain_y',
                    'pdb_aa_y', 'Functional', 'd', 'g',
                    '1/d', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M',
                    'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', *aa, 'pdb_aa_x', 'normalized_score',
                    'total_nan_neigh', 'zr4s_JC'] \
                   + [f for f in huang_pssm_df.columns.tolist()
                      if f.startswith('neighbor_pos_') or f.startswith("n2v")]

    categorical_features = ['glycosylation', 'protein_interaction', 'disorder', 'binding', 'catalysis', 'pdb_aa',
                            'aa_group_5', 'aa_group_HP', 'structure']

    all_features = [f for f in huang_pssm_df.columns.tolist() if f not in not_features]
    numeric_features = [f for f in all_features if f not in categorical_features]
    # if scannet_table:
    #     numeric_features = [f for f in all_features if 'scannet' in f]
    #     categorical_features = []

    logging.debug(f'numeric features = {numeric_features}\n')

    logging.debug(f'categorical features = {categorical_features}')

    all_features_pssm = [*numeric_features, *categorical_features]
    X_train_pssm = huang_pssm_df[all_features_pssm]
    y_train_pssm = huang_pssm_df[aa]
    logging.debug(f'X_train_pssm shape = {X_train_pssm.shape}\n')
    logging.debug(f'X_train_pssm shape = {y_train_pssm.shape}\n')

    ohe = OneHotEncoder(handle_unknown='ignore')
    numeric_transformer_rf = Pipeline(
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


    preprocessor = Pipeline(steps)
    logging.debug(f'preprocessing data\n')
    preprocessor.fit(X_train_pssm)
    X_train_pssm_proc = preprocessor.transform(X_train_pssm)
    logging.debug(f'preprocessing data finished\n')
    logging.debug(f'transformed shape {X_train_pssm_proc.shape}')
    logging.debug(f'new features={preprocessor.named_steps["preprocessor"].transformers_[1][1].named_steps["onehot"].get_feature_names(categorical_features)}')
    model = KerasClassifier(build_fn=NNModel_class,input_dim=[X_train_pssm_proc.shape[1]],num_classes=y_train_pssm.shape[1],vote=True,n_hidden=N_HIDDEN,epochs=10, batch_size=128,verbose=1) #,
    logging.debug(f'ANN chosen\nparam = {model.get_params()}')
    logging.debug('\n')
    logging.debug(f'training')
    model.fit(X_train_pssm_proc, y_train_pssm)
    logging.debug(f'training finished\n')
    logging.debug(f'predicting pssm for training\n')

    PSSM_PRED = model.model.predict(X_train_pssm_proc)

    X_test_proc = preprocessor.transform(test_data[[*numeric_features, *categorical_features]])
    y_hat = model.model.predict(X_test_proc)
    logging.debug(y_hat)

    logging.debug(f'prediction finished\n')
    logging.debug(f'pred shape = {PSSM_PRED}\n')
    for i, AA in enumerate(aa):
        huang_pssm_df[AA + '_pred'] = PSSM_PRED[:,i]

    logging.debug(f'writing predictions to csv file\n')
    huang_pssm_df.to_csv(TRAININIG_PREDICTIONS_OUTPATH,index=False)
    logging.debug(f'saving preprocessor and model\n')
    pickle.dump(preprocessor, open(FITTED_PREPROCESSOR_OUTPATH, 'wb'))
    pickle.dump(model, open(TRAINED_MODEL_OUTPATH, 'wb'))

    logging.debug(f'done\n')
else:
    logging.debug(f'predicting {TEST_DATA_PATH}')
    PREPROCESSOR_PSSM = pickle.load(open(FITTED_PREPROCESSOR_OUTPATH, 'rb'))
    # TRAINED_MODEL_PSSM = load_model(TRAINED_MODEL_OUTPATH)
    TRAINED_MODEL_PSSM = pickle.load(open(TRAINED_MODEL_OUTPATH,'rb'))
    # TRAINED_MODEL_PSSM = load_model("/groups/pupko/natannag/consurf_n2v/huang/trained_pssm_model/trained_ANN_pssm_280722.pkl")
    test_data = pd.read_csv(TEST_DATA_PATH)
    numeric_features = ['2_clique', '3_clique', '4_clique', '5_clique', '6_clique', '7_clique', 'average_neighbor_degree', 'betweenness_centrality', 'clustering_coefficient', 'degree_centrality', 'eigenvector_centrality', 'graphlet1', 'graphlet10', 'graphlet11', 'graphlet12', 'graphlet13', 'graphlet14', 'graphlet15', 'graphlet16', 'graphlet17', 'graphlet18', 'graphlet19', 'graphlet2', 'graphlet20', 'graphlet21', 'graphlet22', 'graphlet23', 'graphlet24', 'graphlet25', 'graphlet26', 'graphlet27', 'graphlet28', 'graphlet29', 'graphlet3', 'graphlet30', 'graphlet31', 'graphlet32', 'graphlet33', 'graphlet34', 'graphlet35', 'graphlet36', 'graphlet37', 'graphlet38', 'graphlet39', 'graphlet4', 'graphlet40', 'graphlet41', 'graphlet42', 'graphlet43', 'graphlet44', 'graphlet45', 'graphlet46', 'graphlet47', 'graphlet48', 'graphlet49', 'graphlet5', 'graphlet50', 'graphlet51', 'graphlet52', 'graphlet53', 'graphlet54', 'graphlet55', 'graphlet56', 'graphlet57', 'graphlet58', 'graphlet59', 'graphlet6', 'graphlet60', 'graphlet61', 'graphlet62', 'graphlet63', 'graphlet64', 'graphlet65', 'graphlet66', 'graphlet67', 'graphlet68', 'graphlet69', 'graphlet7', 'graphlet70', 'graphlet71', 'graphlet72', 'graphlet73', 'graphlet8', 'graphlet9', 'median_rsa_neigh', 'median_wcn_ca_neigh', 'median_wcn_sc_neigh', 'node_degree', 'rsa', 'total_A_neigh', 'total_Aliphatic_neigh', 'total_Aromatic_neigh', 'total_C_neigh', 'total_Charged_neigh', 'total_D_neigh', 'total_Diverse_neigh', 'total_E_neigh', 'total_F_neigh', 'total_G_neigh', 'total_H_neigh', 'total_Hydrophobic_neigh', 'total_I_neigh', 'total_K_neigh', 'total_L_neigh', 'total_M_neigh', 'total_N_neigh', 'total_P_neigh', 'total_Polar_neigh', 'total_Q_neigh', 'total_R_neigh', 'total_S_neigh', 'total_T_neigh', 'total_Tiny_neigh', 'total_V_neigh', 'total_W_neigh', 'total_Y_neigh', 'total_catalytic_neigh', 'total_contact_neigh', 'total_disordered_neigh', 'total_glycosylated_neigh', 'total_interface_neigh', 'total_site_neigh', 'total_ss_B_neigh', 'total_ss_E_neigh', 'total_ss_G_neigh', 'total_ss_H_neigh', 'total_ss_I_neigh', 'total_ss_P_neigh', 'total_ss_S_neigh', 'total_ss_T_neigh', 'total_ss_nan_neigh', 'wcn_ca', 'wcn_sc', 'predicted_score', 'blosum0', 'blosum1', 'blosum2', 'blosum3', 'blosum4', 'blosum5', 'blosum6', 'blosum7', 'blosum8', 'blosum9', 'blosum10', 'blosum11', 'blosum12', 'blosum13', 'blosum14', 'blosum15', 'blosum16', 'blosum17', 'blosum18', 'blosum19', 'aa_freq_ExPASy']
    categorical_features = ['glycosylation', 'protein_interaction', 'disorder', 'binding', 'catalysis', 'pdb_aa', 'aa_group_5', 'aa_group_HP', 'structure']

    X_test_proc = PREPROCESSOR_PSSM.transform(test_data[[*numeric_features, *categorical_features]])
    y_hat = TRAINED_MODEL_PSSM.model.predict(X_test_proc)
    aa = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    for i, AA in enumerate(aa):
        test_data[AA + '_pred'] = y_hat[:, i]
    logging.debug(y_hat)
    test_data.to_csv(RESULTS_DIR+'/'+os.path.split(TEST_DATA_PATH)[-1].split("_")[0]+'_features_and_predictions.csv',index=False)
    test_data.to_csv(RESULTS_DIR+'/'+JOB_TITLE+'_features_and_predictions.csv',index=False)
    with open(RESULTS_DIR+'/'+JOB_TITLE+"_done_pssm_pred.flag",'w') as f:
        f.write('done')