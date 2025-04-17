#

laptop = False
import sys
import os
if laptop:
    path2proteinnet = 'proteinnet/casp12/'
    list_files = ['proteinnet/casp12/validation']
else:
    path2proteinnet_j = '/home/iscb/wolfson/jeromet/EvolutionPrediction/casp12/'
    path2proteinnet_n = '/groups/pupko/natannag/EvolutionPrediction/casp12/'
    path2proteinnet_n2 = '/groups/pupko/natannag/EvolutionPrediction/proteinnet/code'
    if os.path.exists(path2proteinnet_j):
        list_files = [path2proteinnet_j + x for x in ['training_95','validation']]
        path2proteinnet = path2proteinnet_j

    else:
        list_files = [path2proteinnet_n + x for x in ['training_30','validation']]
        sys.path.append(path2proteinnet_n)
        sys.path.append(path2proteinnet_n2)
        path2proteinnet = path2proteinnet_n


# sys.path.append('proteinnet/proteinnet_master/code/')
import text_parser
# from proteinnet.proteinnet_master.code import text_parser
import numpy as np
import pickle





aa = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
      'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V',  'W', 'Y', '-']


def num2seq(sequence):
    return ''.join([aa[x] for x in sequence])







for file in list_files:
    set_name = file.split('/')[-1]
    output_name = path2proteinnet + set_name + '_labels.data'

    list_origins = []
    list_sequences = []
    list_resids = []
    list_PSSMs = []



    with open(file,'r') as file_handle:
        while True:
            try:
                record = text_parser.read_record(file_handle,20)
                try:
                    origin = record['id'].split('#')[-1].split('_')
                    if len(origin) ==3:
                        origin = origin[0].lower() + '_0-' + origin[2]
                    else:
                        origin = origin[0].lower() + '_0-' + origin[1][5:-1].upper()
                    sequence = num2seq(record['primary'])
                    pssm = np.array(record['evolutionary']).T
                    list_origins.append(origin)
                    list_sequences.append(sequence)
                    list_PSSMs.append(pssm)
                except Exception as e:
                    print(record['id'])
                    print(e)
                    continue
            except Exception as e2:
                print(e2)
                break

    env = {
    'list_origins': list_origins,
    'list_sequences':list_sequences,
    'list_resids': list_resids,
    'list_PSSMs' : list_PSSMs
    }
    pickle.dump(env,open(output_name,'wb'),2)






