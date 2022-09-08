import click
@click.command()
@click.argument('dssp_output_dir', type=click.Path(exists=True))
@click.argument('calc_rsa_output_dir', type=click.Path(exists=False))
def main(dssp_output_dir,calc_rsa_output_dir):
    for f in os.listdir(dssp_output_dir):
        input_f = os.path.join(dssp_output_dir,f)
        cmd = f'python calc_rsa.py {input_f} {calc_rsa_output_dir}'
        print(cmd)
        subprocess.check_output(cmd,shell=True)

if __name__ == "__main__":
    import os
    import subprocess
    import argparse
    import csv
    import textwrap
    main()

