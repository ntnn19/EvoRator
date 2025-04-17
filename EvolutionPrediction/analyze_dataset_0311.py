import pandas as pd
import numpy as np
from keras import Sequential
from keras.layers import Dense,BatchNormalization,Activation

import sys

laptop = True
if laptop:
    path2github = '../'
else:
    path2github = '/Users/natan/Documents'
sys.path.append(path2github)

from ScanNet_dev.preprocessing import sequence_utils
def MLP(input_size,nlayers=2,layer_size=32):
    model = Sequential(input_shape=(input_size))
    for k in range(nlayers):
        model.add(Dense(layer_size,name='dense_%s'%k))
        model.add(BatchNormalization())



# table = pd.read_csv('/Users/jerometubiana/Documents/GitHub/EvolutionPrediction/results/model_selection/deep_learning_for_proteinnet/pp_deep_learning_for_proteinnet_n_hidden_30_test_pred.csv')
import datatable as dt
table = pd.read_csv(r"C:\Users\natan\Documents\missense_pathogenecity\source\del_pred_pssm_var_evorator.csv")
# table = dt.fread(r"C:\Users\natan\Documents\missense_pathogenecity\source\del_pred_pssm_var_evorator.csv").to_pandas()
# table = dt.fread(r"C:\Users\natan\Documents\EvolutionPrediction\results\huang\ANN\debug_cwv_ANN_n_hidden_147_test_pred.csv").to_pandas()
# table = dt.fread(r"C:\Users\natan\Documents\EvolutionPrediction\del_pred_pssm_var_evorator.csv").to_pandas()
table[['%s_freq'%a for a in sequence_utils.aa[:20]]] = table[['%s_freq'%a for a in sequence_utils.aa[:20]]].fillna(0)
print(table.columns)
print(table.shape)
print(table['pdb_aa'].unique())
table['pdb_aa'] = table['pdb_aa'].fillna(table['pdb_aa_x'].fillna(table['pdb_aa_y']))
# table['pdb_aa'] = table['pdb_aa'].fillna(table['pdb_aa_x'])
# table['pdb_aa'] = table['pdb_aa'].fillna(table['pdb_aa_y'])
print(table['pdb_aa'].value_counts())
print(table[table['pdb_aa']==''].pdb_position)
print(table.shape)
# table=table[table['pdb_aa']!='']
table=table.dropna(subset=['pdb_aa'])
print(table.shape)
# exit()
original_aa = np.array(table['pdb_aa'])
numerical_original_aa = np.array([sequence_utils.aa.index(x) for x in original_aa])
actual_pssm = np.array(table[['%s_freq'%a for a in sequence_utils.aa[:20]]])
predicted_pssm = np.array(table[['%s_freq_pred'%a for a in sequence_utils.aa[:20]]])

empirical_frequencies_aa = np.array([(numerical_original_aa==k).mean() for k in range(20)])
empirical_frequencies_pssm = actual_pssm.mean(0)
for k in range(20):
    print(sequence_utils.aa[k],empirical_frequencies_aa[k],empirical_frequencies_pssm[k])
accs=[]
ps=[]
for pdb_id in table.pdb_id.unique():
    original_aa = np.array(table[table.pdb_id==pdb_id]['pdb_aa'])
    numerical_original_aa = np.array([sequence_utils.aa.index(x) for x in original_aa])
    actual_pssm = np.array(table[table.pdb_id==pdb_id][['%s_freq' % a for a in sequence_utils.aa[:20]]])
    predicted_pssm = np.array(table[table.pdb_id==pdb_id][['%s_freq_pred' % a for a in sequence_utils.aa[:20]]])

    actual_pssm_ = actual_pssm.copy()
    actual_pssm_[np.arange(len(numerical_original_aa)),numerical_original_aa] = -1
    most_likely_substition = np.argmax(actual_pssm_,axis=-1)

    predicted_pssm_ = predicted_pssm.copy()
    predicted_pssm_[np.arange(len(numerical_original_aa)),numerical_original_aa] = -1
    predicted_most_likely_substition = np.argmax(predicted_pssm_,axis=-1)

    blosum_pssm =  np.array([actual_pssm[numerical_original_aa==k].mean(0) for k in range(20)])
    predicted_pssm_by_blosum = blosum_pssm[numerical_original_aa,:]
    predicted_pssm_by_blosum_ = predicted_pssm_by_blosum.copy()
    predicted_pssm_by_blosum_[np.arange(len(numerical_original_aa)),numerical_original_aa] = -1
    predicted_most_likely_substition_blosum = np.argmax(predicted_pssm_by_blosum_,axis=-1)
    substitution_accuracy = (most_likely_substition == predicted_most_likely_substition).mean()
    substitution_accuracy_blosum = (most_likely_substition == predicted_most_likely_substition_blosum).mean()
    C = np.corrcoef(predicted_pssm.flatten(), actual_pssm.flatten())
    Cblosum = np.corrcoef(predicted_pssm_by_blosum.flatten(), actual_pssm.flatten())
    print(C[0],Cblosum[0])
    print(substitution_accuracy, substitution_accuracy_blosum)
    accs.append((substitution_accuracy, substitution_accuracy_blosum))
    ps.append((C,Cblosum))
##
import matplotlib.pyplot as plt
bins = np.linspace(0, 1, 100)
plt.scatter([acc[0] for acc in accs],[acc[1] for acc in accs])
plt.plot([0,1],[0,1])
# plt.legend(loc='upper right')
plt.show()
##

exit()

accuracy_matrix = np.zeros([20,20])

for i in range(20):
    for j in range(20):
        subset = (numerical_original_aa == i) & (most_likely_substition == j)
        # accuracy_matrix[i,j] = (predicted_most_likely_substition[subset] == j).mean()
        accuracy_matrix[i,j] = (predicted_most_likely_substition_blosum[subset] == j).mean()

import matplotlib.pyplot as plt
import seaborn as sns

plt.figure()
sns.heatmap(accuracy_matrix*100,annot=True,fmt='.0f')
plt.xticks(np.arange(20)+0.5,sequence_utils.aa[:20])
plt.yticks(np.arange(20)+0.5,sequence_utils.aa[:20])
plt.show()

print(most_likely_substition == predicted_most_likely_substition)
print(most_likely_substition == predicted_most_likely_substition_blosum)
substitution_accuracy = (most_likely_substition == predicted_most_likely_substition).mean()
substitution_accuracy_blosum = (most_likely_substition == predicted_most_likely_substition_blosum).mean()

print(substitution_accuracy,substitution_accuracy_blosum)

C = np.corrcoef(predicted_pssm.flatten(),actual_pssm.flatten() )
Cblosum = np.corrcoef(predicted_pssm_by_blosum.flatten(),actual_pssm.flatten())
print()