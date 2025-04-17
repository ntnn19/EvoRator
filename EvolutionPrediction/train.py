import click, sys,os,pickle
laptop = False if os.path.exists("/groups/pupko/natannag/natan_git") else True

def set_num_threads(num_threads=2):
    os.environ["MKL_NUM_THREADS"] = "%s"%num_threads
    os.environ["NUMEXPR_NUM_THREADS"] = "%s"%num_threads
    os.environ["OMP_NUM_THREADS"] = "%s"%num_threads
    os.environ["OPENBLAS_NUM_THREADS"] = "%s"%num_threads
    os.environ["VECLIB_MAXIMUM_THREADS"] = "%s"%num_threads
    os.environ["NUMBA_NUM_THREADS"] = "%s"%num_threads

AA_FREQS = {'A': 8.25 / 100, 'Q': 3.93 / 100, 'L': 9.65 / 100, 'S': 6.63 / 100, 'R': 5.53 / 100,
                'E': 6.72 / 100, 'K': 5.80 / 100, 'T': 5.35 / 100, 'N': 4.05 / 100,
                'G': 7.08 / 100, 'M': 2.41 / 100, 'W': 1.09 / 100, 'D': 5.46 / 100, 'H': 2.27 / 100,
                'F': 3.86 / 100, 'Y': 2.92 / 100, 'C': 1.38 / 100, 'I': 5.91 / 100, 'P': 4.73 / 100,
                'V': 6.86 / 100}


def  calc_conservation_score(PWM):
    # C = log21 +∑logPWM(a)
    PWM = PWM.div(PWM.sum(axis=1), 'index')
    C = np.full((PWM.shape[0],),np.log2(20)) + np.sum(np.log2(PWM.fillna(1e-5)),axis=1)
    print(C)
    return C

# log21 +∑alogPWM(a)

def calc_predicted_pssm_conservation_score(df,AA_FREQS,suffix=''):
    # AA = [i + '_pred' for i in AA_FREQS.keys()]
    AA = [i+suffix for i in AA_FREQS.keys()]


    df[AA] = df[AA].div(df[AA].sum(axis=1), 'index')
    # print(df[AA])
    # print(df[AA].sum(axis=1))
    for aa in AA:
        # df[aa + '_norm_by_freq'] = df[aa] * np.log2(df[aa] / AA_FREQS[aa[0]])
        df[aa + '_norm_by_freq'] = df[aa] * np.log2((df[aa]+1e-5) / AA_FREQS[aa[0]])
    # print(df[AA])
    # print(df[[aa+'_norm_by_freq' for aa in AA]])
    if suffix:
        df[f'{suffix[1:]}_pssm_conservation_score'] = df[[aa + '_norm_by_freq' for aa in AA]].sum(axis=1)
        return
    df[f'pssm_conservation_score'] = df[[aa + '_norm_by_freq' for aa in AA]].sum(axis=1)
    # return df[[aa + '_norm_by_freq' for aa in AA]].sum(axis=1)

def PWM(pred, mut):
    aa_freqs = {'A': 8.25 / 100, 'Q': 3.93 / 100, 'L': 9.65 / 100, 'S': 6.63 / 100, 'R': 5.53 / 100,
                'E': 6.72 / 100, 'K': 5.80 / 100, 'T': 5.35 / 100, 'N': 4.05 / 100,
                'G': 7.08 / 100, 'M': 2.41 / 100, 'W': 1.09 / 100, 'D': 5.46 / 100, 'H': 2.27 / 100,
                'F': 3.86 / 100, 'Y': 2.92 / 100, 'C': 1.38 / 100, 'I': 5.91 / 100, 'P': 4.73 / 100,
                'V': 6.86 / 100}
    return pred / aa_freqs[mut[0]]


def NNModel_class(input_dim=None, num_classes=2, n_hidden=20):
    X_input = Input(input_dim)
    X = Dense(input_dim[0], activation='relu', name='fc2',kernel_regularizer=l2(0.0005))(X_input)
    X = BatchNormalization(name = 'bn0')(X)
    X = Dense(input_dim[0], activation='relu', name='fc3',kernel_regularizer=l2(0.0005))(X)
    X = BatchNormalization(name = 'bn1')(X)
    X = Dense(num_classes, activation='softmax', name='fc_out')(X) # as in preprint
    model = Model(inputs=X_input, outputs=X)
    model.compile('adam', 'kullback_leibler_divergence', metrics=['kullback_leibler_divergence'])  # as in preprint
    return model

def NNModel_reg(input_dim=None, num_classes=2, n_hidden=20):
    X_input = Input(input_dim)
    # X = Dense(input_dim[0]//30, activation='relu', name='fc2', kernel_regularizer=l2(0.0005))(X_input)
    # X = BatchNormalization(name='bn0')(X)
    # X = Dense(input_dim[0]//50, activation='relu', name='fc3', kernel_regularizer=l2(0.0005))(X)
    # X = BatchNormalization(name='bn1')(X)
    X = Dense(num_classes, activation='linear', name='fc_out')(X_input) # as in preprint
    # X = Dense(num_classes, activation='linear', name='fc_out')(X) # as in preprint
    model = Model(inputs=X_input, outputs=X)
    model.compile('adam', 'mean_squared_error', metrics=['mean_squared_error'])
    return model



def permut_protein(protein_df, p):
    permuted_df = protein_df[protein_df.pdb_id == p][
        ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']].astype(
        np.int8).copy().to_numpy()
    np.random.shuffle(permuted_df)
    logging.debug(permuted_df.shape)
    return pd.DataFrame(permuted_df,
                        columns=['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T',
                                 'V', 'W', 'Y']).astype(np.int8)




@click.command()
@click.argument('training_data', type=click.Path(exists=True))
@click.argument('outdir', type=click.Path(exists=None))
@click.option('--parameters', '-p', type=click.Path(exists=True), # you can ignore this parameter for now
              help='Provide a path to parameter file in the following format: <1st parameter> <n-th parameter>\nParameters: ; ANN - n_hidden (size of the hidden layer) ; GBT - n_estimators, max_leaf_nodes, learning rate ; SVC - C, gamma ; LR - l1_ratio')
@click.option('--output-name', '-o', type=str, default='final_output')
@click.option('--exclude-scannet', '-sn', is_flag=True, # you can ignore this parameter for now
              help="Model is trained without the features extracted by scannet")
@click.option('--exclude-evorator', '-er', is_flag=True, # you can ignore this parameter for now
              help="Model is trained without the features extracted by scannet")
@click.option('--train-only', '-to', is_flag=True, help="Model is trained over all data and saved") # you can ignore this parameter for now
@click.option('--model-selection', '-ms', is_flag=True, help="Data are splitted to train, val, and test sets for selecting the best DL architechture")
@click.option('--pssm-prediction', '-pp', is_flag=True, help="PSSM prediction task")
@click.option('--conservation-prediction', '-cp', is_flag=True, help="Conservation scores prediction task") # you can ignore this parameter for now
@click.option('--protein-interaction', '-pi', type=click.Path(exists=True), help="whether to map protein interaction annotations") # you can ignore this parameter for now
def model_selection(training_data, outdir, parameters, output_name, exclude_scannet,exclude_evorator, train_only, model_selection,pssm_prediction,conservation_prediction,protein_interaction):
    # '''
    # This program trains a NN (a multi-layer perceptron) to predict PSSM based on protein structural information.
    # (FEEL FREE TO CHANGE THE ARCHITECTURE OF THE NETWORK using the function <NNClass>; NOTE THE USED LOSS IS <KULLBACK_LEIBLER_DIVERGENCE>. THIS CANNOT BE CHANGED)
    # You can write your own model as long you use the <KULLBACK_LEIBLER_DIVERGENCE> as the loss and <SOFTMAX> activation in the output layer.
    #
    # TRAINING_DATA - csv table in which each row corresponds to an amino acid residue
    #
    # Results are saved at OUTDIR:
    # (1) a table of Pearson correlation coefficients for training, validation , and test datasets for the specified parameters.
    # (2) a table with predictions for test data.
    #
    # CMD = "proteinnet\evorator_mapped_to_consurf\evorator_proteinnet_consurf_100001.csv" results/model_selection --model-selection -pp"
    OUTPATH = os.path.join(outdir, 'deep_learning_for_proteinnet')
    os.makedirs(OUTPATH, exist_ok=True)
    if model_selection:
        train_only= False
    elif train_only:
        TRAINED_MODEL_DIR = os.path.join(OUTPATH, 'trained_models')
        os.makedirs(TRAINED_MODEL_DIR, exist_ok=True)
        if parameters:
            suffix = os.path.split(parameters)[-1].split(".")[0]
        else:
            suffix = 'default'
#        TRAINED_MODEL_OUTPATH = os.path.join(TRAINED_MODEL_DIR, output_name  + '_trained_model_del.pkl')
#        TRAINED_MODEL_OUTPATH = os.path.join(TRAINED_MODEL_DIR, output_name  + '_trained_model_del_2')
        TRAINED_MODEL_OUTPATH = os.path.join(TRAINED_MODEL_DIR, output_name  + '_trained_model_del.h5')
        FITTED_PREPROCESSOR_OUTPATH = os.path.join(TRAINED_MODEL_DIR,
                                                   output_name  + '_fitted_preprocessor_del.pkl')


    # python C:\Users\natan\Documents\missense_pathogenecity\output\huang_aa_variance.csv "C:\Users\natan\Documents\missense_pathogenecity\data\huang_with_evorator_features_and_predictions.csv" -a ANN -p "C:\Users\natan\Documents\missense_pathogenecity\model_selection\param\ANN_param\ANN_n_hidden_23.param"
    # LR C:\Users\natan\Documents\missense_pathogenecity\output\huang_aa_variance.csv "C:\Users\natan\Documents\missense_pathogenecity\data\huang_with_evorator_features_and_predictions.csv" -a LR -p "C:\Users\natan\Documents\missense_pathogenecity\model_selection\param\LR_param\LR_l1_ratio_0.2828282828282829.param



    if parameters:

        with open(parameters, 'r') as f:
            lines = f.read().strip().split(' ')

            n_hidden = int(lines[0])
            parameters_suffix = "n_hidden_" + str(n_hidden)
        # parameters= os.path.split(parameters)[-1].split(".")[0]

    if not parameters:  # default parameter values
        n_hidden = 30
        parameters_suffix = "n_hidden_" + str(n_hidden)

    logging.basicConfig(
        filename=OUTPATH + '/' + f'deep_learning_for_consurf_{output_name}_feature_extraction_and_prediction.log',
        level=logging.DEBUG,
        format='%(asctime)s %(levelname)-8s %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S')
    logging.debug('Starting model selection')
    logging.debug('\n')


    PRED_OUTPATH = os.path.join(OUTPATH, output_name + '_' + 'deep_learning_for_proteinnet' + '_test_pred.csv')
    # TRAINED_MODEL_OUTPATH = os.path.join(OUTPATH,output_name+'_'+parameters+'_trained_model_del.pkl')
    # FITTED_PREPROCESSOR_OUTPATH = os.path.join(OUTPATH,output_name+'_'+parameters+'_fitted_preprocessor_del.pkl')
    CONFIGURATION_AUC_OUTPATH = os.path.join(OUTPATH, output_name + '_performance.csv')

    logging.debug('\n')
    logging.debug(f'Loading data from:\n{training_data}')
    logging.debug('\n')
    # load the data
    #all_data = dt.fread(training_data).to_pandas()
#    flg =True
    flg =False
    if flg:
        logging.debug('read first')
        all_data1 = pd.read_csv(training_data.replace(".csv","_split_aa"))
        logging.debug('drop dup first')
        logging.debug(f'all_data {all_data1.shape}')
        all_data1 = all_data1.drop_duplicates(subset=['pdb_position'])
        logging.debug(f'all_data {all_data1.shape}')
        logging.debug('read second')
        all_data2 = pd.read_csv(training_data.replace(".csv","_split_ab"),header=None)
        logging.debug('read third')
        all_data3 = pd.read_csv(training_data.replace(".csv","_split_ac"),header=None)
        # logging.debug('read fourth')
        # all_data4 = pd.read_csv(training_data.replace(".csv","_split_ad"),header=None)
        # logging.debug('read fifth')
        # all_data5 = pd.read_csv(training_data.replace(".csv","_split_ae"),header=None)
        all_data2.columns = all_data1.columns
        logging.debug('drop dup first')
        logging.debug(f'all_data {all_data2.shape}')
        all_data2 = all_data2.drop_duplicates(subset=['pdb_position'])
        logging.debug(f'all_data {all_data2.shape}')

        all_data3.columns = all_data1.columns
        logging.debug(f'all_data {all_data3.shape}')
        all_data3 = all_data3.drop_duplicates(subset=['pdb_position'])
        logging.debug(f'all_data {all_data3.shape}')
        # all_data4.columns = all_data1.columns
        # all_data5.columns = all_data1.columns
        all_data1 ,all_data2, all_data3 = reduce_memory_usage(all_data1), reduce_memory_usage(all_data2), reduce_memory_usage(all_data3)
        logging.debug('concat')
        # all_data = pd.concat([all_data1,all_data2,all_data3,all_data4,all_data5])
        all_data = pd.concat([all_data1,all_data2,all_data3])



    else:
        # all_data = pd.read_csv(training_data)
        all_data = dt.fread(training_data).to_pandas()
        logging.debug(f'all_data {all_data.shape}')
        all_data = all_data.drop_duplicates(subset=['pdb_position'])
        logging.debug(f'all_data {all_data.shape}')

    logging.debug('\n')
    if protein_interaction:
        ppi_map = pickle.load(open(protein_interaction,'rb'))
        all_data['protein_interaction'] = all_data['pdb_position'].map(ppi_map)
    print(all_data.shape)
    # all_data = all_data.drop_duplicates(subset=['pdb_position'])
    all_data = reduce_memory_usage(all_data)
    matrix = matlist.blosum62


    # ++    +/-  ++  ++  +      +   +     +   ++  +     +    ++  +    ++    ++   ++   +    +    +/-  +/-
    aa = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    # aa = [c + "_norm" for c in aa]
    aa_freqs = {'A': 8.25 / 100, 'Q': 3.93 / 100, 'L': 9.65 / 100, 'S': 6.63 / 100, 'R': 5.53 / 100,
                'E': 6.72 / 100, 'K': 5.80 / 100, 'T': 5.35 / 100, 'N': 4.05 / 100,
                'G': 7.08 / 100, 'M': 2.41 / 100, 'W': 1.09 / 100, 'D': 5.46 / 100, 'H': 2.27 / 100,
                'F': 3.86 / 100, 'Y': 2.92 / 100, 'C': 1.38 / 100, 'I': 5.91 / 100, 'P': 4.73 / 100, 'V': 6.86 / 100}
    aa_pred = [i + '_pred' for i in aa_freqs.keys()]
    all_data['pdb_id'] = all_data['pdb_position'].str.split("_").str[1]
    logging.debug(f'training data dimensions = {all_data.shape}')
    # test_data = test_data.rename(columns={'total_ss_loop_neigh': 'total_ss_nan_neigh'})
    # all_data = all_data.dropna(subset=aa)  # complete case analysis
    # all_data = all_data.rename(columns={'total_ss_loop_neigh': 'total_ss_nan_neigh'})
    # blosum_mat = blosum.blosum62
    # blosum_d = {}
    # for AA1 in aa:
    #     tmp_vec = []
    #     for AA2 in aa:
    #         tmp_vec.append(blosum_mat.get((AA1[0], AA2[0]), blosum_mat.get((AA2[0], AA1[0]))))
    #         blosum_d[AA1[0]] = tmp_vec
    # print(blosum_d)
    # blosum_cols = ['blosum' + str(i) for i in range(len(blosum_d['A']))]
    # all_data[blosum_cols] = all_data['pdb_aa'].map(blosum_d).apply(pd.Series)
    # all_data['aa_freq_ExPASy'] = all_data['pdb_aa'].map(aa_freqs)
    # for mut in aa:
    #     all_data[mut] = all_data[mut].apply(PWM, mut=mut)
    # print(blosum_d)
    # exit()

    all_data[[AA[0] for AA in aa]] = all_data[[AA[0] for AA in aa]].fillna(1e-5)  # complete case analysis
    all_data['split_by_CATH']  = all_data['split_by_CATH'].str[:4]+"_"+all_data['chain']
    print(all_data.split_by_CATH.value_counts(normalize=True).sort_values(ascending=False))
    # calc conservation score log20 +∑logPWM(a)
    all_data['pssm_conservation_score'] = calc_predicted_pssm_conservation_score(all_data,aa_freqs)

    # all_data['pssm_conservation_score'] = calc_conservation_score(all_data[aa])
    print(all_data[['pssm_conservation_score','evorator_conservation_score','consurf_conservation_score']].corr())
    all_data_clean = all_data.dropna(subset=[*[AA[0] for AA in aa],'split_by_CATH','consurf_conservation_score'])  # complete case analysis with respect to class



    if model_selection:
        p_list = all_data_clean['split_by_CATH'].unique()

        logging.debug(f'Number of proteins = {len(p_list)}')
        # split to train, val, and test
        random.seed(65)
        train_p_list,val_p_list, test_p_list,not_in_cath = train_test_split.get_splitting_indices(list(p_list))
        all_data_clean= all_data_clean[~all_data_clean['split_by_CATH'].isin(not_in_cath)]
        logging.debug(f'Complete case analysis dimensions= {all_data_clean.shape}')
        logging.debug(f'Train proteins={train_p_list},len={len(train_p_list)}')
        logging.debug(f'Val proteins={val_p_list},len={len(val_p_list)}')
        logging.debug(f'Test proteins={test_p_list},len={len(test_p_list)}')
        logging.debug('\n')
        logging.debug(
            f'Sanity check independent data\n\nExpected set size = 201\nActual set size = {len(list(set([*train_p_list, *test_p_list, *val_p_list])))}')
        logging.debug('\n')

    all_data_clean = all_data_clean.sample(frac=1).reset_index(drop=True)
    not_features = ['index', 'pdb_position', 'pdb_id', 'chain', 'pdb', 'pdb_x', 'pdb_y', 'chain_x', 'chain_y',
                    'pssm_conservation_score', 'consurf_conservation_score', 'confidence_interval', '_merge',
                    'split_by_CATH','total_interface_neigh',
                    'pdb_aa_y', 'Functional', 'd', 'g',
                    '1/d', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M',
                    'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', *aa, 'pdb_aa_x',
                    'total_nan_neigh', 'zr4s_JC'] \
                   + [f for f in all_data_clean.columns.tolist()
                      if f.startswith('neighbor_pos_') or f.startswith("n2v") or f.endswith('_norm_by_freq')]


    print('Adding ScanNet features')
    # scannet_features = scannet_features[['pdb_position',*scannet_features.iloc[:,96+128:].columns.tolist()]]
    scannet_features_names = [f for f in all_data_clean.columns.tolist() if 'scannet_feature' in f]
    print(scannet_features_names)
    print(f'Resulting data size= {all_data_clean.shape}')

    categorical_features = ['glycosylation', 'disorder', 'binding', 'catalysis', 'pdb_aa',
                            'aa_group_5', 'aa_group_HP', 'structure']
    if protein_interaction:
        categorical_features.append('protein_interaction')

    all_features = [f for f in all_data_clean.columns.tolist() if f not in not_features]
    numeric_features = [f for f in all_features if f not in categorical_features]
    features_that_should_be_numeric = np.array(numeric_features)[((all_data_clean[numeric_features].dtypes) == np.dtype('O'))]  # remove later, debugging
    print('numeric_features',all_data_clean[numeric_features])
    for c in numeric_features:
        try:
            all_data_clean[c] = all_data_clean[c].str.decode('utf8')
        except Exception as e:
            print('ERROR')
            print(e)

    print('features_that_should_be_numeric',all_data_clean[numeric_features])
    print('features_that_should_be_numeric',all_data_clean[features_that_should_be_numeric])


    logging.debug(f'numeric features = {numeric_features}\n')

    logging.debug(f'categorical features = {categorical_features}')

    # clean_train = all_data_clean[all_data_clean.pdb_id.isin(train_p_list[:2])]
    if model_selection:
        print('train',all_data_clean[all_data_clean.split_by_CATH.isin(p_list[train_p_list].tolist())]['split_by_CATH'].unique().shape)
        print('val',all_data_clean[all_data_clean.split_by_CATH.isin(p_list[val_p_list].tolist())]['split_by_CATH'].unique().shape)
        print('test',all_data_clean[all_data_clean.split_by_CATH.isin(p_list[test_p_list].tolist())]['split_by_CATH'].unique().shape)
        clean_train = all_data_clean[all_data_clean.split_by_CATH.isin(p_list[train_p_list].tolist())]
        clean_val = all_data_clean[all_data_clean.split_by_CATH.isin(p_list[val_p_list].tolist())]
        clean_test = all_data_clean[all_data_clean.split_by_CATH.isin(p_list[test_p_list].tolist())]

        X_train = clean_train[[*numeric_features, *categorical_features]]
        X_val = clean_val[[*numeric_features, *categorical_features]]
        X_test = clean_test[[*numeric_features, *categorical_features]]

        # clean_train[numeric_features]=clean_train[numeric_features].astype(float)
        if pssm_prediction:
            y_train = clean_train[[AA[0] for AA in aa]].astype(np.float16)
            y_val = clean_val[[AA[0] for AA in aa]].astype(np.float16)
            y_test = clean_test[[AA[0] for AA in aa]].astype(np.float16)
            conservation_prediction = False
        elif conservation_prediction:
            y_train = clean_train['consurf_conservation_score'].astype(np.float16)
            y_val = clean_val['consurf_conservation_score'].astype(np.float16)
            y_test = clean_test['consurf_conservation_score'].astype(np.float16)
            pssm_prediction = False

        logging.debug('\n')
        logging.debug(f'X_train={X_train.shape},y_train={y_train.shape}')
        logging.debug(f'X_val={X_val.shape},y_val={y_val.shape}')
        logging.debug(f'X_test={X_test.shape},y_test={y_test.shape}')

    elif train_only:
        logging.debug('training over all data')
        X_train = all_data_clean[[*numeric_features, *categorical_features]].dropna(how='all', axis=1)
        if pssm_prediction:
            y_train = all_data_clean[[AA[0] for AA in aa]].astype(np.float16)
            conservation_prediction = False
        elif conservation_prediction:
            y_train = all_data_clean['consurf_conservation_score'].astype(np.float16)
            pssm_prediction = False

    ohe = OneHotEncoder(handle_unknown='ignore')
    numeric_transformer_rf = Pipeline(
        steps=[('imputer', SimpleImputer(strategy='median')), ('scaler', StandardScaler()),
               ('variance_filter', VarianceThreshold())])
    categorical_transformer = Pipeline(
        steps=[('imputer', SimpleImputer(strategy='most_frequent')), ('onehot', ohe),
               ('variance_filter', VarianceThreshold())])
    exclude_clique  = False
    if exclude_clique:
        numeric_features = [f for f in numeric_features if not f.endswith('clique')]

    if exclude_scannet:
        print('scannet excluded')
        numeric_features = [f for f in numeric_features if f not in scannet_features_names]
        column_transformer = ColumnTransformer(transformers=[('num', numeric_transformer_rf, numeric_features),('cat', categorical_transformer, categorical_features)])
    elif exclude_evorator:
        print('evorator excluded')
        numeric_features =  scannet_features_names
        column_transformer = ColumnTransformer(transformers=[('num', numeric_transformer_rf, numeric_features),('cat', 'drop', categorical_features)])
        categorical_features =[]

    else:
        column_transformer = ColumnTransformer(transformers=[('num', numeric_transformer_rf, numeric_features),('cat', categorical_transformer, categorical_features)])
    steps = [('preprocessor', column_transformer)]
    print(numeric_features)
    logging.debug(f'numeric {numeric_features}')
    print(categorical_features)
    logging.debug(f'categorical {categorical_features}')
    print(categorical_features)
    # if not os.path.exists(FITTED_PREPROCESSOR_OUTPATH):


    preprocessor = Pipeline(steps)
    preprocessor.fit(X_train)
    # pickle.dump(preprocessor, open(FITTED_PREPROCESSOR_OUTPATH, 'wb'))
    # else:
    #     logging.debug('\n')
    #     logging.debug(f'Found fitted preprocessor at {FITTED_PREPROCESSOR_OUTPATH}')
    #     logging.debug('\n')
    #     preprocessor = pickle.load(open(FITTED_PREPROCESSOR_OUTPATH,'rb'))

    logging.debug(f'Imputing, standardizing, one-hot encoding, and constant variance filtering of features')
    X_train_proc = preprocessor.transform(X_train)
    if model_selection:
        X_val_proc = preprocessor.transform(X_val)
        X_test_proc = preprocessor.transform(X_test)
        logging.debug(f'X_test_proc={X_test_proc.shape},y_test={y_test.shape}')
        logging.debug(f'X_val_proc={X_val_proc.shape},y_val={y_val.shape}')

    logging.debug(f'X_train_proc={X_train_proc.shape},y_train={y_train.shape}')
    try:
        logging.debug(
        f'new features={preprocessor.named_steps["preprocessor"].transformers_[1][1].named_steps["onehot"].get_feature_names(categorical_features)}')
    except Exception as e:
        print(e)

    #
    # if os.path.exists(TRAINED_MODEL_OUTPATH):

    # logging.debug(f'Found trained algorithm at {TRAINED_MODEL_OUTPATH}')
    # logging.debug('\n')
    # model = pickle.load(open(TRAINED_MODEL_OUTPATH,'rb'))

    # model = LogisticRegression(max_iter=5000)
    logging.debug(f'Training with params {parameters}')
    logging.debug('\n')

        # model = NNModel_class([X_train_proc.shape[1]], num_classes=y_train.shape[1],n_hidden=n_hidden)
    if pssm_prediction:
        if model_selection:
            ES = tf.keras.callbacks.EarlyStopping(monitor='val_loss', patience=10)
            model = KerasClassifier(build_fn=NNModel_class, input_dim=[X_train_proc.shape[1]], num_classes=y_train.shape[1],
                                    n_hidden=n_hidden, epochs=50, batch_size=128, verbose=1,validation_data=(X_val_proc, y_val),callbacks=[ES])  # ,
        elif train_only:
            logging.debug('Training over all data')
            SM = tf.keras.callbacks.ModelCheckpoint(filepath=TRAINED_MODEL_OUTPATH,
                                                 save_weights_only=True,
                                                 verbose=1)

            model = KerasClassifier(build_fn=NNModel_class, input_dim=[X_train_proc.shape[1]],
                                    num_classes=y_train.shape[1],
                                    n_hidden=n_hidden, epochs=10, batch_size=128, verbose=1,callbacks=[SM])# ,
        logging.debug(f'ANN chosen\nparam = {model.get_params()}')
        logging.debug('\n')

    elif conservation_prediction:
        ES = tf.keras.callbacks.EarlyStopping(monitor='val_loss', patience=5)
        model = KerasRegressor(build_fn=NNModel_reg, input_dim=[X_train_proc.shape[1]], num_classes=1,
                                n_hidden=n_hidden, epochs=50, batch_size=128, verbose=1,validation_data=(X_val_proc, y_val),callbacks=[ES])  # ,
        logging.debug(f'ANN chosen\nparam = {model.get_params()}')
        logging.debug('\n')

    model.fit(X_train_proc, y_train)
    logging.debug(f'Training finished')
    if train_only:

#        model.model.save(TRAINED_MODEL_OUTPATH)
#         model.model.save(TRAINED_MODEL_OUTPATH)
#         model.model.save_weights(TRAINED_MODEL_OUTPATH)
        pickle.dump(preprocessor, open(FITTED_PREPROCESSOR_OUTPATH, 'wb'))

        exit()


    logging.debug(f'Training finished')
    # Predict training and validation , stack them, refit preprocessor and algorithm and predict for test
    # clean_train_new = all_data_clean[all_data_clean.pdb_id.isin([*train_p_list, *val_p_list][:2])]
    # clean_train_new = all_data_clean[all_data_clean.pdb_id.isin([*train_p_list, *val_p_list])]
    # X_train_new = clean_train_new[[*numeric_features, *categorical_features]]
    # y_train_new = clean_train_new[aa].astype(np.float16)

    pred_training_set = model.model.predict(X_train_proc)
    pred_val_set = model.model.predict(X_val_proc)
    pred_test_set = model.model.predict(X_test_proc)

    if pssm_prediction:
        train_score_to_assert = stats.pearsonr(y_train.values.ravel() ,pred_training_set.ravel())
        val_score_to_assert = stats.pearsonr(y_val.values.ravel() ,pred_val_set.ravel())
        test_score_to_assert = stats.pearsonr(y_test.values.ravel() ,pred_test_set.ravel())
        train_score_spearman = stats.spearmanr(y_train.values.ravel() ,pred_training_set.ravel())
        val_score_spearman = stats.spearmanr(y_val.values.ravel() ,pred_val_set.ravel())
        test_score_spearman = stats.spearmanr(y_test.values.ravel() ,pred_test_set.ravel())
        logging.debug('train',train_score_to_assert)
        logging.debug('val',val_score_to_assert)
        logging.debug('test',test_score_to_assert)
        clean_test = pd.concat([clean_test,pd.DataFrame(pred_test_set,columns=[AA[0] + '_predicted' for AA in aa])],axis=1)
        clean_train = pd.concat([clean_train,pd.DataFrame(pred_training_set,columns=[AA[0] + '_predicted' for AA in aa])],axis=1)
        clean_val = pd.concat([clean_val,pd.DataFrame(pred_val_set,columns=[AA[0] + '_predicted' for AA in aa])],axis=1)
        train_score_after_concat = stats.pearsonr(clean_train[[*[AA[0] + '_predicted' for AA in aa],*aa]][aa].dropna().values.ravel(),clean_train[[*[AA[0] + '_predicted' for AA in aa],*aa]][[AA[0] + '_predicted' for AA in aa]].dropna().values.ravel())
        logging.debug('train_after_concat', train_score_after_concat)
        val_score_after_concat = stats.pearsonr(clean_val[[*[AA[0] + '_predicted' for AA in aa],*aa]][aa].dropna().values.ravel(),clean_val[[*[AA[0] + '_predicted' for AA in aa],*aa]][[AA[0] + '_predicted' for AA in aa]].dropna().values.ravel())
        logging.debug('val_after_concat', val_score_after_concat)
        test_score_after_concat = stats.pearsonr(clean_test[[*[AA[0] + '_predicted' for AA in aa],*aa]][aa].dropna().values.ravel(),clean_test[[*[AA[0] + '_predicted' for AA in aa],*aa]][[AA[0] + '_predicted' for AA in aa]].dropna().values.ravel())
        logging.debug('test_after_concat', test_score_after_concat)
        assert train_score_to_assert ==  train_score_after_concat, f'expected pearson r {train_score_to_assert}, got {train_score_after_concat}'
        assert val_score_to_assert ==  val_score_after_concat, f'expected pearson r {val_score_to_assert}, got {val_score_after_concat}'
        assert test_score_to_assert ==  test_score_after_concat, f'expected pearson r {test_score_to_assert}, got {test_score_after_concat}'
        # clean_test[[AA[0] + '_predicted' for AA in aa]] = pred_test_set
        # clean_train[[AA[0] + '_predicted' for AA in aa]] = pred_training_set
        # clean_val[[AA[0] + '_predicted' for AA in aa]] = pred_val_set

    if conservation_prediction:
        clean_test['evorator2_conservation_score'] = pred_test_set
        clean_train['evorator2_conservation_score'] = pred_training_set
        clean_val['evorator2_conservation_score'] = pred_val_set

    logging.debug('\n')
    logging.debug(np.asarray(pred_training_set).shape)
    logging.debug(np.asarray(pred_test_set).shape)
    logging.debug(np.asarray(pred_val_set).shape)
    if pssm_prediction:
        pass
        # clean_train['pssm_conservation_score'] = calc_conservation_score(clean_train[aa])
        # clean_train['predicted_pssm_conservation_score'] = calc_conservation_score(clean_train[aa_pred])
        #
        # clean_val['pssm_conservation_score']= calc_conservation_score(clean_val[aa])
        # clean_val['predicted_pssm_conservation_score']= calc_conservation_score(clean_val[aa_pred])
        #
        # clean_test['pssm_conservation_score']= calc_conservation_score(clean_test[aa])
        # clean_test['predicted_pssm_conservation_score']= calc_conservation_score(clean_test[aa_pred])

        # calc_predicted_pssm_conservation_score(clean_train,AA_FREQS)
        # calc_predicted_pssm_conservation_score(clean_val,AA_FREQS)
        # calc_predicted_pssm_conservation_score(clean_test,AA_FREQS)
        #
        # calc_predicted_pssm_conservation_score(clean_train, AA_FREQS, suffix='_predicted')
        # calc_predicted_pssm_conservation_score(clean_val, AA_FREQS, suffix='_predicted')
        # calc_predicted_pssm_conservation_score(clean_test, AA_FREQS, suffix='_predicted')
        #
        # print(clean_train[['pssm_conservation_score','predicted_pssm_conservation_score','evorator_conservation_score','consurf_conservation_score']].corr())
        # print(clean_test[['pssm_conservation_score','predicted_pssm_conservation_score','evorator_conservation_score','consurf_conservation_score']].corr())
        # print(clean_val[['pssm_conservation_score','predicted_pssm_conservation_score','evorator_conservation_score','consurf_conservation_score']].corr())
        #
        metric= 'n_largest_accuracy' # average across proteins
        # metric1= 'n_largest_accuracy' # average across proteins
        # metric2= 'grouoed_n_largest_accuracy' # average across proteins
        score_train = train_score_to_assert,train_score_spearman
        score_val = val_score_to_assert,val_score_spearman
        score_test = test_score_to_assert,test_score_spearman
        # score_train = clean_train[["pssm_conservation_score","predicted_pssm_conservation_score"]].corr().iloc[-1,0]
        # score_val = clean_val[["pssm_conservation_score","predicted_pssm_conservation_score"]].corr().iloc[-1,0]
        # score_test = clean_test[["pssm_conservation_score","predicted_pssm_conservation_score"]].corr().iloc[-1,0]
        # sns.jointplot(x=clean_test[[*[AA[0] + '_predicted' for AA in aa],*aa]][aa].values.ravel(),y=clean_test[[*[AA[0] + '_predicted' for AA in aa],*aa]][[AA[0] + '_predicted' for AA in aa]].values.ravel(), kind="scatter", color="#4CB391")
        #sns.jointplot(x=y_test.values.ravel() ,y=pred_test_set.ravel(), kind="hist", color="#4CB391")
        #plt.show()

    if conservation_prediction:
        # average across proteins
        score_train = r2_score(y_train,pred_training_set)
        score_val = r2_score(y_val,pred_val_set)
        score_test = r2_score(y_test,pred_test_set)
        metric = 'R2'

    logging.debug('\n')
    logging.debug(f'Train {metric} = {score_train}')
    logging.debug('\n')
    logging.debug('\n')
    logging.debug(f'Validation {metric} = {score_val}')
    logging.debug('\n')
    logging.debug('\n')
    logging.debug(f'Test {metric} = {score_test}')
    logging.debug('\n')
    results_df = pd.DataFrame([[score_train, score_val, score_test]],
                              columns=[f'{metric}_train', f'{metric}_val', f'{metric}_test'])
    results_df.to_csv(CONFIGURATION_AUC_OUTPATH, index=False)

    clean_test.to_csv(PRED_OUTPATH, index=False)
    logging.debug('Completed!')



if __name__ == '__main__':
    import time
    import numpy as np
    import matplotlib.pyplot as plt
    import pandas as pd
    import random
    from sklearn.impute import SimpleImputer
    import os
    import logging
    import pickle
    from sklearn.compose import ColumnTransformer

    from sklearn.metrics import  r2_score
    from sklearn.impute import SimpleImputer
    from sklearn.compose import ColumnTransformer
    from sklearn.pipeline import  Pipeline
    from sklearn.feature_selection import VarianceThreshold
    from sklearn.preprocessing import  OneHotEncoder, StandardScaler

    from tensorflow.python.keras.layers import Input, Dense, Activation, Flatten, Conv2D, BatchNormalization
    from tensorflow.python.keras.wrappers.scikit_learn import  KerasClassifier, KerasRegressor
    from tensorflow.python.keras.models import Model
    from tensorflow.python.keras.regularizers import l2
    from Bio.SubsMat import MatrixInfo as matlist
    import pickle
    import random
    import train_test_split
    import seaborn as sns
    import datatable as dt
    import tensorflow as tf
    from reduce_memory_usage import reduce_memory_usage
    from scipy import stats
    # set_num_threads(4)
    model_selection()

