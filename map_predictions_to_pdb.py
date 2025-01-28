import pandas as pd
from biopandas.pdb import PandasPdb


def map_scores_to_pdb(score_file, pdb_file, output_pdb):
    """
    Maps scores from a table to a PDB file by updating the B-factor column.

    Parameters:
    - score_file (str): Path to the CSV file containing residue scores.
    - pdb_file (str): Path to the PDB file.
    - output_pdb (str): Path to save the modified PDB file.
    """

    # Load score table
    scores_df = pd.read_csv(score_file)
    print(scores_df.head())

    # Load PDB file
    pdb = PandasPdb().read_pdb(pdb_file)

    # Extract ATOM records
    atoms_df = pdb.df['ATOM']

    # Parse residue IDs from the score table (extract chain and residue number)
    scores_df[['Chain', 'Residue']] = scores_df['pdb_position'].str.extract(r'([A-Z])(\d+)', expand=True)
    scores_df['Residue'] = scores_df['Residue'].astype(int)  # Convert residue number to integer

    # Filter for Cα atoms only
    atoms_df = atoms_df.merge(scores_df, left_on=['chain_id', 'residue_number'], right_on=['Chain', 'Residue'],
                              how='left')

    # Replace B-factor column with scores, filling NaNs with 0
    atoms_df['b_factor'] = atoms_df['Score'].fillna(0)
    print(atoms_df.b_factor.describe())
    # Save the modified PDB file
    pdb.df['ATOM'] = atoms_df
    pdb.to_pdb(path=output_pdb, records=['ATOM', 'HETATM'], gz=False)

    print(f"Modified PDB file saved to: {output_pdb}")
# Example usage
# map_scores_to_pdb('scores.csv', 'input.pdb', 'output.pdb')
