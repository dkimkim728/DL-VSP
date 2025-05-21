import pandas as pd
import numpy as np
import os

# -------------------------------
# Log Fold Change Calculation
# -------------------------------
# Description:
#   Computes log₂ fold change between treatment and control 
#   expression values for the gene of interest.
# -------------------------------

# --- File paths (same output structure as prior steps) ---
input_dir = "1_data_processing/outputs/"
output_file = os.path.join(input_dir, "final_training_data.csv")

treatment_file = os.path.join(input_dir, "final_treatment_avg.csv")
control_file = os.path.join(input_dir, "final_control_avg.csv")

# --- Load datasets ---
treatment_df = pd.read_csv(treatment_file)
control_df = pd.read_csv(control_file).rename(columns={"gene_expression": "control_expression"})

# --- Merge on cell line and calculate log₂ fold change ---
merged = pd.merge(treatment_df, control_df[["cell_mfc_name", "control_expression"]],
                  on="cell_mfc_name", how="left")

merged["LFC"] = np.log2(merged["gene_expression"]) - np.log2(merged["control_expression"])

# --- Drop invalid or incomplete rows ---
final_df = merged.dropna(subset=["LFC", "canonical_smiles"])

# --- Save output ---
final_df.to_csv(output_file, index=False)
print(f"✅ Final training data saved to: {output_file}")