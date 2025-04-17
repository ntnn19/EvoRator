import numpy as np
import os
import pandas as pd
#import datatable as dt
import sys
import logging
from reduce_memory_usage import reduce_memory_usage

TABLE_DIR = sys.argv[1]
DATE = sys.argv[2] # DDMMYY



logging.basicConfig(
        filename=TABLE_DIR + '/' + f'concat_evorator_scannet_consurf_pssm_{DATE}.log',
        level=logging.DEBUG,
        format='%(asctime)s %(levelname)-8s %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S')


dfs = []
i=0
for f in os.listdir(TABLE_DIR):
    try:
        if i%100==0:
        # if i==100:
            print(i,":",f)
            logging.debug(f"{i}:{f}")
            # break

        #df = dt.fread(os.path.join(TABLE_DIR,f)).to_pandas()
        df = pd.read_csv(os.path.join(TABLE_DIR,f))

        # df['split_by_CATH'] = f.split("_")[0].lower()+df['chain']
        df['split_by_CATH'] = f.split("_")[0].lower()[:4]+"_"+df['chain']
        df = reduce_memory_usage(df,verbose=False)
        dfs.append(df)
        i+=1

    except Exception as e:
        print(e)
        logging.debug(e)


final_df = pd.concat(dfs,axis=0)

final_df.to_csv(os.path.join(TABLE_DIR,"evorator_scannet_pssm_consurf_{DATE}.csv"),index=False)
#
