#!/usr/bin/python

'''
This script parses an input PDB file and returns relative solvent accessibility
values (RSA), secondary structure information, and raw DSSP output.

Author: Benjamin R. Jack
'''


import click

# ASA normalization constants were taken from:
# M. Z. Tien, A. G. Meyer, D. K. Sydykova, S. J. Spielman, C. O. Wilke (2013).
# Maximum allowed solvent accessibilities of residues in proteins. PLOS ONE
# 8:e80635.
RES_MAX_ACC = {'A': 129.0, 'R': 274.0, 'N': 195.0, 'D': 193.0, \
               'C': 167.0, 'Q': 225.0, 'E': 223.0, 'G': 104.0, \
               'H': 224.0, 'I': 197.0, 'L': 201.0, 'K': 236.0, \
               'M': 224.0, 'F': 240.0, 'P': 159.0, 'S': 155.0, \
               'T': 172.0, 'W': 285.0, 'Y': 263.0, 'V': 174.0}


def parse_dssp_line(line):
    '''
    Extract values from a single line of DSSP output and calculate RSA.
    '''
    solvent_acc = line[35:39].strip()  # record SA value for given AA
    amino_acid = line[13].strip()  # retrieve amino acid
    residue = line[6:10].strip()
    chain = line[11].strip()
    secondary_structure = line[16].strip() # secondary structure
    if amino_acid.islower():
        # if the AA is a lowercase letter, then it's a cysteine
        amino_acid = "C"
    if amino_acid in RES_MAX_ACC:
        max_acc = RES_MAX_ACC[amino_acid]  # Find max SA for residue
        rsa = float(solvent_acc) / max_acc # Normalize SA
    else:
        rsa = None
    return {'pdb_position': residue,
            'pdb_aa': amino_acid,
            'chain': chain,
            'rsa': rsa,
            'structure': secondary_structure}

def parse_dssp(raw_dssp_output):
    '''
    Parse a DSSP output file and return a dictionary of amino acids, PDB residue
    number, chain, and RSA.
    '''
    with open(raw_dssp_output, 'r') as dssp_file:
        lines = dssp_file.readlines()
        output = []  # list of dictionaries for output
    for line in lines[28:]:
        # Skip first 28 lines of header info
        output_line = parse_dssp_line(line)
        if output_line['rsa'] is not None:
            # Skip lines with no RSA value
            output.append(output_line)

    return output


@click.command()
@click.argument('dssp_output_file', type=click.Path(exists=True))
@click.argument('calc_rsa_output_dir', type=click.Path(exists=False))
def main(dssp_output_file,calc_rsa_output_dir):
    '''
    Run mkdssp on input PDB and parse output into CSV with RSA values.
    '''


    output_rsa = os.path.join(calc_rsa_output_dir,os.path.split(dssp_output_file)[-1] + '.rsa.csv')
    asa_file = dssp_output_file


    # DSSP succeeded, write output to CSV
    output_dict = parse_dssp(asa_file)
    with open(output_rsa, 'w') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=['pdb_position',
                                                     'chain', 'structure',
                                                     'pdb_aa', 'rsa'])
        writer.writeheader()
        writer.writerows(output_dict)

if __name__ == "__main__":
    import os
    import subprocess
    import argparse
    import csv
    import textwrap
    main()
