import CATH_utils
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

def get_splitting_indices(list_pdb_ids = ['1a3x_A','1a3x_B','1ntk_A', '2kho_A','13pk_A']):

    list_cath_identifiers = [CATH_utils.cath_database[pdb] for pdb in list_pdb_ids if pdb in list(CATH_utils.cath_database.keys())]
    not_in_cath = [pdb for pdb in list_pdb_ids if pdb not in list(CATH_utils.cath_database.keys())]
    list_cath_identifiers = [[u[0] for u in v] for v in list_cath_identifiers]
    nexamples = len(list_pdb_ids)-len(not_in_cath)

    shares_cath = csr_matrix( (nexamples,nexamples),dtype=np.bool)
    for u in range(nexamples):
        for v in range(nexamples):
            if any([x in list_cath_identifiers[u] for x in list_cath_identifiers[v]]):
                shares_cath[u,v] = 1
    n_components, component_labels = connected_components(shares_cath,directed=False)

    train_indices = []
    validation_indices = []
    test_indices = []

    for component in range(n_components):
        subset = list(np.nonzero(component_labels == component)[0])
        if component % 10 == 0:
            test_indices += subset
        elif component % 15 == 1:
            validation_indices += subset
        else:
            train_indices += subset

    train_indices = np.array(train_indices)
    validation_indices = np.array(validation_indices)
    test_indices = np.array(test_indices)

    return train_indices, validation_indices, test_indices, not_in_cath



