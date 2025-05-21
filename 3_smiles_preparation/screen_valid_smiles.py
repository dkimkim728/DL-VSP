# screen_smiles_validity.py

from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd
import os

# === INPUT ===
input_file = "input_smiles.csv"     # must contain a 'SMILES' column
output_valid_csv = "validated_smiles.csv"
mol_output_dir = "mol2d_preview"
os.makedirs(mol_output_dir, exist_ok=True)

# === LOAD AND VALIDATE SMILES ===
df = pd.read_csv(input_file)
valid_smiles = []
valid_mols = []

for smiles in df["SMILES"]:
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None and mol.GetNumAtoms() > 3:  # complexity filter
        valid_smiles.append(smiles)
        valid_mols.append(mol)

# === SAVE VALID SMILES TO CSV ===
pd.DataFrame({"SMILES": valid_smiles}).to_csv(output_valid_csv, index=False)

# === OPTIONAL: EXPORT FIRST 100 MOLECULES TO MOL FILES ===
for i, mol in enumerate(valid_mols[:100]):
    AllChem.Compute2DCoords(mol)
    Chem.MolToMolFile(mol, os.path.join(mol_output_dir, f"mol_{i}.mol"))

print(f"✅ {len(valid_smiles)} valid SMILES saved to {output_valid_csv}")
print(f"🧪 2D structures exported for top {min(100, len(valid_mols))} molecules in '{mol_output_dir}/'")