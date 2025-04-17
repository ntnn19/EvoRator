import os
mode="portable"
cwd=os.getcwd()
print("CWD=",cwd)
if mode=="portable":
    library_folder = cwd+'ScanNet_dev/'
    structures_folder = library_folder + 'PDB/'# Where pdb/mmCIF structures files are stored.
    MSA_folder = library_folder+'MSA/' # Where multiple sequence alignments are stored.
    predictions_folder = library_folder + 'predictions/' # Output folder.
    model_folder = library_folder + 'models/' # Where the networks as stored as pairs of files (.h5,.data).
    pipeline_folder = library_folder + 'pipelines/' # Where preprocessed datasets are stored.
    initial_values_folder = model_folder + 'initial_values/' # Where initial values of the parameters for the gaussian kernels and residue-residue graph edges are stored.
    homology_folder = library_folder + 'baselines/homology/' # Where files are stored for homology baseline.
    visualization_folder = library_folder + 'visualizations/'
    path2hhblits = library_folder # Path to hhblits binary. Not required if using ScanNet_noMSA networks.
    path2sequence_database = None # Path to sequence database Not required if using ScanNet_noMSA networks.
    path_to_dssp = '' # Path to dssp binary. Only for reproducing baseline performance.
    path_to_msms = '' # Path to msms binary. Only for reproducing baseline performance.
    path_to_multiprot = None # Path to multiprot executable. Only relevant for homology baseline
    path_to_cdhit = ''
    path_to_mafft = 'mafft'
