"""
mol2_to_frcmod.py

Generates FRCMOD parameter files from GAFF-formatted MOL2 ligand files using Antechamber's parmchk2.
"""

import os
import subprocess

# === Directory Configuration ===
input_dir = "5_md_simulation/inputs/gaff_mol2"   # Input: GAFF MOL2 files
output_dir = "5_md_simulation/inputs/frcmod"     # Output: FRCMOD files

# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

# === Process All .mol2 Files ===
mol2_files = [f for f in os.listdir(input_dir) if f.endswith('.mol2')]

for mol2_file in mol2_files:
    input_path = os.path.join(input_dir, mol2_file)
    output_path = os.path.join(output_dir, mol2_file.replace(".mol2", ".frcmod"))

    print(f"🔁 Processing: {mol2_file}")
    try:
        result = subprocess.run(
            ["parmchk2", "-i", input_path, "-f", "mol2", "-o", output_path],
            check=True, capture_output=True, text=True
        )

        if result.stderr:
            print(f"⚠️ stderr: {result.stderr}")
        print(f"✅ Generated: {os.path.basename(output_path)}")

    except subprocess.CalledProcessError as e:
        print(f"❌ Failed: {mol2_file} → {e}")
    except Exception as e:
        print(f"❌ Unexpected error: {mol2_file} → {e}")