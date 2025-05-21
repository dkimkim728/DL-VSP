"""
convert_pdb_to_mol2.py

Converts a ligand PDB file to a GAFF-compatible MOL2 file using RDKit, Open Babel, and Antechamber.
Input:  ligand_name.pdb  in  5_md_simulation/inputs/pdb/
Output: GAFF-formatted MOL2 in  5_md_simulation/inputs/gaff_mol2/
"""

import os
import subprocess
from rdkit import Chem
from rdkit.Chem import AllChem

# === Ligand Identifier ===
ligand_id = "ligand001"  # <-- Change this to your ligand file prefix (e.g., ligand001 for ligand001.pdb)

# === Directory Setup ===
root_dir = "5_md_simulation/inputs"
pdb_dir = os.path.join(root_dir, "pdb")
mol2_fixed_dir = os.path.join(root_dir, "mol2_fixed")
gaff_mol2_dir = os.path.join(root_dir, "gaff_mol2")

# === Ensure Output Folders Exist ===
os.makedirs(mol2_fixed_dir, exist_ok=True)
os.makedirs(gaff_mol2_dir, exist_ok=True)

# === Define File Paths ===
pdb_file = os.path.join(pdb_dir, f"{ligand_id}.pdb")
mol_file = os.path.join(mol2_fixed_dir, f"{ligand_id}_rdkit.mol")
mol2_file = os.path.join(mol2_fixed_dir, f"{ligand_id}_rdkit.mol2")
gaff_mol2_file = os.path.join(gaff_mol2_dir, f"{ligand_id}_gaff.mol2")

print(f"\n🔁 Processing ligand: {ligand_id}")

# === Step 1: Load and Prepare Molecule with RDKit ===
mol = Chem.MolFromPDBFile(pdb_file, sanitize=True, removeHs=False)
if mol is None:
    raise ValueError(f"❌ Failed to load ligand from: {pdb_file}")

mol = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol, AllChem.ETKDG())
AllChem.UFFOptimizeMolecule(mol)
AllChem.ComputeGasteigerCharges(mol)

Chem.MolToMolFile(mol, mol_file)
print(f"✅ RDKit MOL saved: {mol_file}")

# === Step 2: Convert MOL → MOL2 (Open Babel) ===
subprocess.run(["obabel", mol_file, "-O", mol2_file], check=True)
print(f"✅ Converted to MOL2: {mol2_file}")

# === Step 3: Generate GAFF MOL2 (Antechamber) ===
subprocess.run([
    "antechamber",
    "-i", mol2_file, "-fi", "mol2",
    "-o", gaff_mol2_file, "-fo", "mol2",
    "-at", "gaff", "-c", "bcc",
    "-nc", "0", "-m", "1", "-s", "2"
], check=True)
print(f"✅ GAFF MOL2 saved: {gaff_mol2_file}")