import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, Lipinski

# Mute RDKit parsing warnings
RDLogger.DisableLog('rdApp.*')

# Initialize collections
valid_mols = []
exclusion_log = []
table_rows = []
seen_canonical = set()

def truncate_smiles(smiles, max_len=20):
    """Truncates SMILES strings and adds dots if they exceed max_len."""
    return f"{smiles[:max_len]}..." if len(smiles) > max_len else smiles

def lipinski_violations(mol):
    """Calculates the number of Lipinski's Rule of 5 violations for a molecule."""
    logp = Descriptors.MolLogP(mol)
    h_donors = Lipinski.NumHDonors(mol)
    h_acceptors = Lipinski.NumHAcceptors(mol)
    
    return logp, h_donors, h_acceptors

def process_smiles(smiles_list):
    for raw_s in smiles_list:
    # Sanitization and parsing
        mol = Chem.MolFromSmiles(raw_s)

        if mol is None:
            exclusion_log.append({
                "Identifier": truncate_smiles(raw_s, 30),
                "Pipeline_Stage": "SMILES Sanitization",
                "Rejection_Reason": "Invalid SMILES Syntax / Unparsable Graph"
            })
            continue

        # Deduplication
        canonical_s = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True) # Allows enantiomers to co-exist
        if canonical_s in seen_canonical:
            exclusion_log.append({
                "Identifier": truncate_smiles(raw_s, 30),
                "Pipeline_Stage": "Deduplication",
                "Rejection_Reason": f"Duplicate structure normalized to: {truncate_smiles(canonical_s)}"
            })
            continue
            
        seen_canonical.add(canonical_s)

        lipinski_data = lipinski_violations(mol)
        logp, h_donors, h_acceptors = lipinski_data

        exp_lookup = dict(zip(df['smiles'], df['exp']))
        exp_val = exp_lookup.get(raw_s, None)

        valid_mols.append(mol)
            
        table_rows.append({
            "SMILES": canonical_s,
            "LogP": round(logp, 2),
            "H-Donors": h_donors,
            "H-Acceptors": h_acceptors,
            "exp": exp_val,
            "TPSA": Descriptors.TPSA(mol),
            "Rotatable_Bonds": Descriptors.NumRotatableBonds(mol),
            "FpDensity": Descriptors.FpDensityMorgan2(mol),
            "FracCSP3": Descriptors.FractionCSP3(mol),
            "MolWt": Descriptors.MolWt(mol)
        })

def main(smiles_list):
    process_smiles(smiles_list)
    
    df_clean = pd.DataFrame(table_rows)
    df_exclusion_log = pd.DataFrame(exclusion_log)

    print(f"\nPipeline Processing Complete.")
    print(f"Total Source Rows: {len(smiles_list)} | Valid Compounds: {len(df_clean)} | Total Rejected Elements: {len(df_exclusion_log)}\n")

    if not df_clean.empty:
        output_csv_path = DATA_PATH / "clean_lipophilicity.csv"
        df_clean.to_csv(output_csv_path, index=False)
        print(f"Cleaned dataset saved to {output_csv_path}")
    else:
        print("No valid compounds found.")

    print("\nEXCLUSION LOG")
    if not df_exclusion_log.empty:
        print(df_exclusion_log.to_string(index=False))
    else:
        print("No compounds were excluded.")
    
    return df_clean
    
df_features_final = main(df.smiles)
