import pandas as pd
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

# load data
input_file = "../intermediate/active_compounds.tsv"
df = pd.read_csv(input_file, sep='\t')

def process_smiles(smiles):
    if pd.isna(smiles):
        return None, None
    
    try:
        mol_input = Chem.MolFromSmiles(smiles)
        if mol_input is None:
            return None, None
        canonical_input = Chem.MolToSmiles(mol_input, isomericSmiles=True)
        
        # sanitize
        sanitized_output = rdMolStandardize.StandardizeSmiles(smiles)
        
    except Exception as e:
        return None, None

    return canonical_input, sanitized_output

# apply
results = df['SMILES'].apply(process_smiles)

# results
df['SMILES_canonical_input'] = [res[0] for res in results]
df['SMILES_sanitized'] = [res[1] for res in results]

# compare
df['changed'] = df['SMILES_canonical_input'] != df['SMILES_sanitized']

num_changed = df['changed'].sum()
total_mols = len(df)

print(f"molecules processed: {total_mols}")
print(f"SMILES changed: {num_changed}")

# format
output_cols = ['AID', 'CID', 'SMILES', 'SMILES_sanitized', 'InChI', 'InChIKey']
df[output_cols].to_csv("../intermediate/active_compounds_sanitized.tsv", sep='\t', index=False)

print("Results saved to '../intermediate/active_compounds_sanitized.tsv'")
