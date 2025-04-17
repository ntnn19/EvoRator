#
import sys

laptop = False
if laptop:
    path2github = '../'
else:
    path2github = '/home/iscb/wolfson/jeromet/'


import numpy as np
import pickle



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
    path2github_j = '/home/iscb/wolfson/jeromet/'
    path2github_n = '/groups/pupko/natannag/'
    if os.path.exists(path2proteinnet_j):
        list_files = [path2proteinnet_j + x for x in ['training_95','validation']]
        path2github=path2github_j

    else:
        list_files = [path2proteinnet_n + x for x in ['training_30','validation']]
        sys.path.append(path2proteinnet_n)
        sys.path.append(path2proteinnet_n2)
        path2github=path2github_n
sys.path.append(path2github)


from ScanNet_dev.utilities import io_utils,sequence_alignment_utils
from ScanNet_dev.preprocessing import PDBio, PDB_processing

for file in list_files:
    env = pickle.load(open(file+'_labels.data','rb'),encoding='latin1')
    list_origins_ = env['list_origins']
    list_sequences_ = env['list_sequences']
    list_resids_ = env['list_resids']
    list_PSSMs_ = env['list_PSSMs']
    print(list_origins_[0])

    list_origins = []
    list_sequences = []
    list_resids = []
    list_PSSMs = []

    L = len(list_origins_)
    for l in range(L):
        origin = list_origins_[l]
        sequence = list_sequences_[l]
        PSSM = list_PSSMs_[l]
        RESID_IDX = list_resids_[l]
        try:
            file, chain_ids = PDBio.getPDB(origin,biounit=False)
            chains = PDBio.load_chains(file=file,chain_ids=chain_ids)[1]
            sequence_from_pdb = PDB_processing.process_chain(chains)[0]
            residue_indices = PDB_processing.get_PDB_indices(chains,return_model=True,return_chain=True)
            mapping_seq2pdb, mapping_pdb2seq, similarity = sequence_alignment_utils.align_and_map(sequence,sequence_from_pdb, return_similarity=True,local=True)
            correct_mapping = (similarity > 0.97)
            if correct_mapping:
                sequence = ''.join([sequence[x] for x in mapping_seq2pdb])
                PSSM = PSSM[mapping_seq2pdb]
                resids = residue_indices[mapping_pdb2seq]
            else:
                sequence = sequence_from_pdb
                PSSM = np.nan * np.ones([len(sequence),20])
                resids = residue_indices
        except Exception as e:
            print(e)
            continue
        list_origins.append(origin)
        list_sequences.append(sequence)
        list_PSSMs.append(PSSM.astype(np.float32))
        list_resids.append(resids)

    env['list_origins'] = list_origins
    env['list_sequences'] = list_sequences
    env['list_resids'] = list_resids
    env['list_PSSMs'] = list_PSSMs
    io_utils.save_pickle(env,file+'_labels.data')





