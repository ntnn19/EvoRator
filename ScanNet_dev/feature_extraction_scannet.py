import os
import re
import sys
import pandas as pd
import numpy as np
#sys.path.append("../ScanNet_dev")
#from ScanNet_dev import predict_features
import predict_features
# sys.path.append("../missense_pathogenicity/casp12")
sys.path.append("/groups/pupko/natannag/EvolutionPrediction/proteinnet/code")
sys.path.append("/groups/pupko/natannag/casp12")
PDB_PATH = sys.argv[1]
OUTPUT_DIR= sys.argv[2]


for f in os.listdir(PDB_PATH):
    if f.endswith('pdb'):
        list_chains = [
                       os.path.join(PDB_PATH,f).replace("\\","/")
                       ]
        print(list_chains)
        PDB_ID = os.path.split(list_chains[0])[-1].split(".")[0]
    
        OUTPUT_NAME = PDB_ID+"_scannet_features.csv"
        OUTPATH = os.path.join(OUTPUT_DIR,OUTPUT_NAME)
    
        layers = ['SCAN_filters_atom_aggregated_activity', # dim 64  aggergated activity of atomic neighborhood
            'all_embedded_attributes_aa', # dim 96 (small neighborhood)
            'SCAN_filter_activity_aa', # dim 128 (larger neighborhood ~13A)
            'SCAN_filters_aa_embedded_1', # dim 32 (low-dimensional projection from SCAN_filter_activity_aa)
        ]
    
        list_dictionary_features = predict_features.predict_features(list_chains,
                                                                     model='ScanNet_PPI_noMSA', # PPBS model without evolution.
                                                                     layer=layers, # layers[0], # AA-scale spatio-chemical filters
                                                                     output_format='dictionary',
                                                                     permissive=False)
    
        # consurf_df['pdb_position'] = consurf_df['chain'].astype(str) +  consurf_df['pdb_position_consurf'].astype(str) + "_" + job_title
        mapping_evorator_key = [k[1]+str(k[-1])+"_"+ PDB_ID for k in list_dictionary_features[0].keys()]
        # print(mapping_evorator_key)
        # exit()
        # print('Residue ID' ,'Features 1-5')
        # print(list_dictionary_features[0].keys())
        # print(len(list_dictionary_features[0].keys()))
        # print(list_dictionary_features[0][(0, 'A', 7)][0].shape)
        # print(list_dictionary_features[0][(0, 'A', 7)][1].shape)
        # print(list_dictionary_features[0][(0, 'A', 7)][2].shape)
        # print(np.concatenate(list_dictionary_features[0][(0, 'A', 7)]))
        # print(np.concatenate(list_dictionary_features[0][(0, 'A', 7)]).shape)
        features_per_res = [np.concatenate(list_dictionary_features[0][k]) for k in list_dictionary_features[0]]
        FINAL_DF = pd.DataFrame(features_per_res,index=mapping_evorator_key,columns=['scannet_feature'+str(i) for i in range(len(features_per_res[0]))]).reset_index()
        FINAL_DF = FINAL_DF.rename(columns={'index':'pdb_position'})
        FINAL_DF.to_csv(OUTPATH,index=False)
    # print(list_dictionary_features.keys())
    # for key ,item in list(list_dictionary_features[0].items())[:10]:
    #     # print(key ,'%.2f,%.2f,%.2f,%.2f,%.2f ' %(item[0] ,item[1] ,item[2] ,item[3] ,item[4] ) )
    #     print(key ,'%.2f,%.2f,%.2f'%(item[0] ,item[1] ,item[2] ) )
