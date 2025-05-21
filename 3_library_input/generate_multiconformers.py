import os
import gc
import subprocess
import pandas as pd
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers
from rdkit.Chem.MolStandardize.rdMolStandardize import TautomerEnumerator

# === Configuration ===
input_csv_path = "./input/ligands.csv"
expected_smiles_column = "SMILES"
output_csv_path = "./output/variants_summary.csv"
store_path = "./output/conformer_files"
rmsd_threshold = 0.5
num_conformers = 10

os.makedirs(store_path, exist_ok=True)

# === Variant Generators ===
def generate_stereoisomers(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return []
    isomers = list(EnumerateStereoisomers(mol))
    return [Chem.MolToSmiles(iso) for iso in isomers]

def generate_tautomers(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return smiles
    tautomer_enumerator = TautomerEnumerator()
    tautomer = tautomer_enumerator.Canonicalize(mol)
    return Chem.MolToSmiles(tautomer)

# === Conformer Generator ===
def generate_conformers(smiles, compound_dir, compound_name):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False
    mol = Chem.AddHs(mol)
    AllChem.EmbedMultipleConfs(mol, numConfs=num_conformers * 2)
    unique_confs = []

    for conf in mol.GetConformers():
        if not unique_confs:
            unique_confs.append(conf)
            continue
        is_unique = all(
            rdMolAlign.GetBestRMS(mol, mol, prbId=conf.GetId(), refId=prev.GetId()) > rmsd_threshold
            for prev in unique_confs
        )
        if is_unique:
            unique_confs.append(conf)
        if len(unique_confs) >= num_conformers:
            break

    if not unique_confs:
        return False

    for i, conf in enumerate(unique_confs):
        pdb_path = os.path.join(compound_dir, f"{compound_name}_conf{i}.pdb")
        pdbqt_path = os.path.join(compound_dir, f"{compound_name}_conf{i}.pdbqt")
        Chem.MolToPDBFile(mol, pdb_path, confId=conf.GetId())
        cmd = f"obabel {pdb_path} -O {pdbqt_path} --gen3D --partialcharge eem --ff GAFF"
        try:
            subprocess.run(cmd, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except subprocess.CalledProcessError:
            continue
        os.remove(pdb_path)
    return True

# === Main Loop ===
def process_smiles():
    output_data = []
    df = pd.read_csv(input_csv_path)
    if expected_smiles_column not in df.columns:
        raise ValueError("SMILES column not found!")

    for idx, row in tqdm(df.iterrows(), total=len(df), desc="Processing Compounds"):
        smiles = row[expected_smiles_column]
        compound_name = f"cmpd_{idx}"
        compound_dir = os.path.join(store_path, compound_name)
        os.makedirs(compound_dir, exist_ok=True)

        stereoisomers = generate_stereoisomers(smiles)
        tautomers = list({generate_tautomers(s) for s in stereoisomers})
        variants = list(set(stereoisomers + tautomers))

        for v_idx, variant in enumerate(variants):
            variant_name = f"{compound_name}_v{v_idx}"
            success = generate_conformers(variant, compound_dir, variant_name)
            if success:
                output_data.append([variant_name, variant])

        gc.collect()

    pd.DataFrame(output_data, columns=["Variant", "SMILES"]).to_csv(output_csv_path, index=False)
    print(f"✅ Finished. Summary saved to {output_csv_path}")

if __name__ == "__main__":
    process_smiles()