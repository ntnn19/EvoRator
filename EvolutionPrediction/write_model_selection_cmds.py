OUTPUT = 'train.cmds'
sizes = [1e5,2e5,3e5,4e5,5e5,6e5,7e5,8e5,9e5,1e6]
tasks = ['pp','cp']
excluded_features = ['-sn','-er',None]
# i=0
configurations = []
for n in sizes:
    for task in tasks:
        for excluded_feature in excluded_features:
            configurations.append([int(n),task,excluded_feature])

with open(OUTPUT,'w') as f:
    for i, config in enumerate(configurations):
        if None in config:
            n,task = config[0], config[1]
            cmd = f'module load python/python-anaconda3.7-itaym;source activate /groups/pupko/natannag/conda/envs/NatanEnv; python /bioseq/evorator/EvolutionPrediction/train.py /bioseq/evorator/curr_consurf_db_outputs/evorator_scannet_pssm_consurf_{n}.csv /bioseq/evorator/ml_results/model_selection -ms -{task} -o debug_{n}_{task}_ppi -pi /bioseq/evorator/consurf_db_random/ppi_map.data\tdebug_{n}_{task}_{i}\n'
        else:
            n,task, excluded_feature = config[0], config[1], config[2]
            cmd = f'module load python/python-anaconda3.7-itaym;source activate /groups/pupko/natannag/conda/envs/NatanEnv; python /bioseq/evorator/EvolutionPrediction/train.py /bioseq/evorator/curr_consurf_db_outputs/evorator_scannet_pssm_consurf_{n}.csv /bioseq/evorator/ml_results/model_selection -ms -{task} -o debug_{n}_{task}_{excluded_feature[1:]}_exclued_ppi {excluded_feature} -pi /bioseq/evorator/consurf_db_random/ppi_map.data\tdebug_{n}_{task}_{i}\n'

        f.write(cmd)