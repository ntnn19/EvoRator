import click

@click.command()
@click.argument('output', type=click.Path(exists=False))
def main(output):
    '''module load python/python-anaconda3.7-itaym;python evorator_2_predict_aa.py /groups/pupko/natannag/consurf_n2v/huang/proteins_features_final/1YSCA/1YSC_features.csv /groups/pupko/natannag/consurf_n2v/huang/proteins_features_final/1YSCA/consurf_summary.txt /groups/pupko/natannag/consurf_n2v/huang/proteins_features_final/1YSCA/ 1YSCA'''

    SCRIPT_PATH= '/groups/pupko/natannag/consurf_n2v/huang/evorator_2_predict_aa_150322_2.py'
    FEATURES_PATH = f'/groups/pupko/natannag/missense_pathogenicity/adress_db/'
    CONSURF_OUTPUT_PATTERN = '/groups/pupko/natannag/missense_pathogenicity/adress_consurf/'
    CMDS_OUTPUT_PATTERN = 'evorator_mapping_residue_variety_cmds.cmds'
    OUTPUT = os.path.join(output,CMDS_OUTPUT_PATTERN)
    fo = open(OUTPUT,'w')
    flg=False
    i=0
    consurf_output_2_path_D = {}
    for f in os.listdir(CONSURF_OUTPUT_PATTERN):
        pdb_id = f.split(".")
        consurf_query = pdb_id[0]
        consurf_output_2_path_D[consurf_query] = [os.path.join(FEATURES_PATH,consurf_query[:-1],consurf_query[:-1]+'_features_and_predictions.csv'),
            os.path.join(CONSURF_OUTPUT_PATTERN,f),
            os.path.join(FEATURES_PATH,consurf_query[:-1])]

    for k,v in consurf_output_2_path_D.items():

        # cmd = f'module load python/python-anaconda3.7-itaym;python {SCRIPT_PATH} {features_path} {consurf_output_path} {subdir} {pdb_id}\tevo{i}\n'
        cmd = f'/groups/pupko/natannag/conda/envs/NatanEnv/bin/python {SCRIPT_PATH} {v[0]} {v[1]} {v[2]} {k[:-1]}\tevo{i}\n'
        fo.write(cmd)
        print(i)
        i += 1


    fo.close()



    # for i,f in enumerate(os.listdir(EVORATOR_PATH)):
    #     if f.endswith(EVORATOR_OUTPUT_PATTERN):
    #         evorator_path_tmp = os.path.join(EVORATOR_PATH,f)
    #         cmd= f'module load python/python-anaconda3.7-itaym;python {SCRIPT_PATH} {missense_path} {evorator_path_tmp} {results_dir} --data={data} --gene-name={pdb_id}\tevo{i}\n'
    #         fo.write(cmd)

if __name__ == '__main__':
    import itertools
    import time
    import numpy as np
    import pandas as pd
    import os
    import logging
    import pickle
    import re
    main()
