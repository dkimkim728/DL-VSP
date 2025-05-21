"""
run_docking_batch.py

Batch docking script using AutoDock Vina.
Processes multiple ligand conformers stored in subdirectories and saves docking results per compound.
"""

import os
import subprocess

# === Configurable Input Paths ===
config_path = "docking_inputs/conf.txt"                    # Docking box and parameters
receptor_path = "docking_inputs/receptor.pdbqt"            # Target protein (converted PDBQT)
compounds_dir = "docking_inputs/conformer_pdbqt/"          # Folder with compound subfolders
output_dir = "4_docking/outputs/"                          # Save docking results here

# === Create Output Directory If Not Exists ===
os.makedirs(output_dir, exist_ok=True)

# === Batch Docking Process ===
for compound in os.listdir(compounds_dir):
    compound_path = os.path.join(compounds_dir, compound)

    if os.path.isdir(compound_path):
        compound_output_path = os.path.join(output_dir, compound)
        os.makedirs(compound_output_path, exist_ok=True)

        for fname in os.listdir(compound_path):
            if fname.endswith(".pdbqt"):
                ligand_path = os.path.join(compound_path, fname)
                output_path = os.path.join(compound_output_path, fname)

                vina_cmd = [
                    "vina",
                    "--config", config_path,
                    "--ligand", ligand_path,
                    "--receptor", receptor_path,
                    "--out", output_path
                ]

                print(f"⏳ Docking: {ligand_path}")
                try:
                    subprocess.run(vina_cmd, check=True)
                    print(f"✅ Success: {output_path}")
                except subprocess.CalledProcessError as e:
                    print(f"❌ Failed: {ligand_path}\n{e}")