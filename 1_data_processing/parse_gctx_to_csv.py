import os
import pandas as pd
from cmapPy.pandasGEXpress import parse

# -------------------------------
# Generalized GCTX Parsing Script
# -------------------------------
# Description:
#   Extract expression values for a specified gene ID from a GCTX file
#   and merge them with selected metadata fields.
# -------------------------------

# --- Placeholder file paths (replace before use) ---
treatment_gctx = "data/your_treatment_data.gctx"
control_gctx = "data/your_control_data.gctx"
metadata_path = "data/instinfo_metadata.txt"

output_dir = "1_data_processing/outputs/"
os.makedirs(output_dir, exist_ok=True)

treatment_csv = os.path.join(output_dir, "filtered_treatment.csv")
control_csv = os.path.join(output_dir, "filtered_control.csv")

# --- Metadata fields to include ---
selected_columns = [
    "sample_id",
    "pert_iname",
    "cell_mfc_name",
    "pert_dose", "pert_dose_unit",
    "pert_itime", "pert_time_unit",
    "pert_type"
]

# --- Function: Extract expression values for one gene ID ---
def extract_expression_table(gctx_file, metadata_df, output_file, gene_id):
    gct = parse.parse(gctx_file)
    expression = gct.data_df                      # Genes × Samples
    row_metadata = gct.row_metadata_df
    col_metadata = gct.col_metadata_df

    # Confirm gene ID is present
    if gene_id not in row_metadata["id"].values:
        raise ValueError(f"Gene ID {gene_id} not found in GCTX file: {gctx_file}")

    gene_expr = expression.loc[row_metadata["id"] == gene_id].squeeze()
    gene_expr.name = "gene_expression"

    # Match sample IDs and merge with metadata
    metadata_filtered = metadata_df[
        metadata_df["sample_id"].isin(gene_expr.index)
    ].copy()
    metadata_filtered = metadata_filtered[selected_columns].set_index("sample_id")
    metadata_filtered["gene_expression"] = gene_expr[metadata_filtered.index]

    metadata_filtered.reset_index().to_csv(output_file, index=False)
    print(f"✅ Saved filtered data to: {output_file}")

# --- Main script block ---
if __name__ == "__main__":
    # Replace this with your gene ID of interest
    gene_id_of_interest = 5300  # Example: PIN1 = 5300

    # Load metadata (TSV format from CLUE or LINCS)
    metadata_df = pd.read_csv(metadata_path, sep="\t", low_memory=False)

    # Run extraction for both treatment and control GCTX files
    extract_expression_table(treatment_gctx, metadata_df, treatment_csv, gene_id_of_interest)
    extract_expression_table(control_gctx, metadata_df, control_csv, gene_id_of_interest)