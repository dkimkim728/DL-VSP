import pandas as pd
import os

# -------------------------------
# Expression Table Preprocessing
# -------------------------------
# Description:
#   - Averages replicate instances by compound and cell line
#   - Appends SMILES information using a compound info file
# -------------------------------

# --- Paths (matching parse_gctx_to_csv.py) ---
input_dir = "1_data_processing/outputs/"
output_dir = "1_data_processing/outputs/"
os.makedirs(output_dir, exist_ok=True)

# Input files
treatment_csv = os.path.join(input_dir, "filtered_treatment.csv")
control_csv = os.path.join(input_dir, "filtered_control.csv")
compound_info_file = "data/compoundinfo_metadata.txt"

# Output files
output_treatment_avg = os.path.join(output_dir, "final_treatment_avg.csv")
output_control_avg = os.path.join(output_dir, "final_control_avg.csv")

# --- Function: average replicate expression by compound & cell line ---
def average_replicates(df):
    avg_df = df.groupby(['pert_iname', 'cell_mfc_name'], as_index=False)['gene_expression'].mean()
    return avg_df

def main():
    # Load expression tables
    print("🔹 Loading filtered treatment and control datasets...")
    treat_df = pd.read_csv(treatment_csv)
    control_df = pd.read_csv(control_csv)

    print("🔹 Averaging replicate instances...")
    treat_avg = average_replicates(treat_df)
    control_avg = average_replicates(control_df)

    # Load compound info with SMILES
    print("🔹 Loading compound metadata and merging SMILES...")
    compound_info = pd.read_csv(compound_info_file, sep="\t", low_memory=False)
    compound_smiles = compound_info[['pert_iname', 'canonical_smiles']].drop_duplicates()

    # Merge SMILES into both tables
    treat_final = pd.merge(treat_avg, compound_smiles, on='pert_iname', how='left')
    control_final = pd.merge(control_avg, compound_smiles, on='pert_iname', how='left')

    # Save to output files
    treat_final.to_csv(output_treatment_avg, index=False)
    control_final.to_csv(output_control_avg, index=False)

    print(f"✅ Final treatment file saved to: {output_treatment_avg}")
    print(f"✅ Final control file saved to: {output_control_avg}")

if __name__ == "__main__":
    main()