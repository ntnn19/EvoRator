import logging

import click

def is_aa(res):
    protein_letters_3 = [k.upper() for k in IUPACData.protein_letters_3to1.keys()]
    if res in protein_letters_3:
        return True
    return False

def get_disorder(parser, pdb_file, pdb_chain):
    l=[]
    structure_id = pdb_file.split('/')[-1].split('.')[0].upper()
    structure = parser.get_structure(structure_id, pdb_file)
    model = structure[0]
    chain = model[pdb_chain]
#    l.append([pdb_chain + str(i.get_full_id()[3][1]) + "_" + structure_id + pdb_chain for i in chain if i.is_disordered() == 1 and is_aa(i.get_resname())])
    l.append([pdb_chain + str(i.get_full_id()[3][1]) + "_" + structure_id.split("_")[0] for i in chain if i.is_disordered() == 1 and is_aa(i.get_resname())])
    l= [res for sublist in l for res in sublist]
    return l

def get_binding_sites(pdb_file):
    l=[]
    structure_id = pdb_file.split('/')[-1].split('.')[0].upper()
    pdb_file = open(pdb_file,'r')
    for i, line in enumerate(pdb_file):
        if line.startswith("SITE     "):
            tmp= [line.strip().split()[4:][i:i+3] for i in range(0,len(line.strip().split()[4:]),3) if is_aa(line.strip().split()[4:][i:i+3][0])]
#            tmp= ["".join(i[1:]) + "_" + structure_id + i[1]  for i in tmp]
            tmp= ["".join(i[1:]) + "_" + structure_id.split("_")[0]  for i in tmp]
            l.append(tmp)
    l = [res for sublist in l for res in sublist]
    return l


@click.command()
@click.argument('pdb_file', type=click.Path(exists=True))
@click.argument('pdb_chain', type=str)
@click.argument('output_dir', type=click.Path(exists=True))
@click.option('--job-title',type=str,default='',show_default=True,help='Insert job title')
def main(pdb_file, pdb_chain,output_dir,job_title):

    parser = PDBParser(PERMISSIVE=1)
    disordered_residues =get_disorder(parser, pdb_file, pdb_chain)
    binding_residues= get_binding_sites(pdb_file)
    disordered_out = open(os.path.join(output_dir,job_title+".disorder.txt"),'w')
    binding_out = open(os.path.join(output_dir,job_title+".binding.txt"),'w')
    disordered_out.write("\n".join(disordered_residues))
    binding_out.write("\n".join(binding_residues))
    disordered_out.close()
    binding_out.close()



if __name__ == '__main__':
    import time
    import pandas as pd
    import numpy as np
    import os
    import requests
    import bs4
    from bs4 import BeautifulSoup
    import urllib
    import subprocess
    import re
    import networkx
    from collections import Counter
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
    from Bio.PDB import PDBParser
    from Bio.PDB.DSSP import DSSP
    from Bio.SeqUtils import IUPACData
    main()
