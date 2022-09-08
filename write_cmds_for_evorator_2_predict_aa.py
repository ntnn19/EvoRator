import click

@click.command()
@click.argument('output', type=click.Path(exists=False))
def main(output):
    '''module load python/python-anaconda3.7-itaym;python evorator_2_predict_aa.py /groups/pupko/natannag/consurf_n2v/huang/proteins_features_final/1YSCA/1YSC_features.csv /groups/pupko/natannag/consurf_n2v/huang/proteins_features_final/1YSCA/consurf_summary.txt /groups/pupko/natannag/consurf_n2v/huang/proteins_features_final/1YSCA/ 1YSCA'''

    SCRIPT_PATH= '/groups/pupko/natannag/consurf_n2v/huang/evorator_2_predict_aa.py'
    FEATURES_PATH = f'/groups/pupko/natannag/consurf_n2v/huang/proteins_features_final'
    CMDS_OUTPUT_PATTERN = 'evorator_mapping_residue_variety_cmds.cmds'
    OUTPUT = os.path.join(output,CMDS_OUTPUT_PATTERN)
    fo = open(OUTPUT,'w')
    flg=False
    i=0
    for subdir, dirs, files in os.walk(FEATURES_PATH):
        for f in files:
            if f.endswith('_features.csv'):
                if os.path.exists(os.path.join(subdir, 'consurf_summary.txt')):
                    print(os.path.join(subdir, f))
                    print(os.path.join(subdir, 'consurf_summary.txt'))
                    features_path = os.path.join(subdir, f)
                    consurf_output_path = os.path.join(subdir, 'consurf_summary.txt')
                    pdb_id = f.split("_")[0]
                    cmd = f'module load python/python-anaconda3.7-itaym;python {SCRIPT_PATH} {features_path} {consurf_output_path} {subdir} {pdb_id}\tevo{i}\n'
                    fo.write(cmd)
                    print(i)
                    i+=1
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
