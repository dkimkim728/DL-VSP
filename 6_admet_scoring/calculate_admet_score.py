import pandas as pd

# === Paths ===
input_csv = '~/results.csv'
output_csv = '~/admet_score.csv'

# === Properties and weights from Guan et al. ===
property_names = [
    'AMES', 'Caco2_Wang', 'Carcinogens_Lagunin', 'CYP1A2_Veith', 'CYP2C19_Veith',
    'CYP2C9_Veith', 'CYP2C9_Substrate_CarbonMangels', 'CYP2D6_Veith',
    'CYP2D6_Substrate_CarbonMangels', 'CYP3A4_Veith', 'CYP3A4_Substrate_CarbonMangels',
    'hERG', 'HIA_Hou', 'Pgp_Broccatelli'
]

weights = [
    0.6021, 0.3074, 0.6857, 0.2379, 0.2685, 0.2966, 0.2633,
    0.2829, 0.2296, 0.3421, 0.2691, 0.4059, 0.7548, 0.3735
]

# === Load input CSV ===
df = pd.read_csv(input_csv)

# === Check if required columns exist ===
required_cols = ['SMILES'] + property_names
if not all(col in df.columns for col in required_cols):
    raise ValueError(f"Input CSV must contain the following columns: {required_cols}")

# === Define substitution function ===
def substitute_values(row):
    substituted = []
    for prop in property_names:
        val = row[prop]
        if prop in ['AMES', 'Carcinogens_Lagunin', 'CYP1A2_Veith', 'CYP2C19_Veith',
                    'CYP2C9_Veith', 'CYP2D6_Veith', 'CYP3A4_Veith', 'hERG', 'HIA_Hou', 'Pgp_Broccatelli']:
            substituted.append(1 if 0 <= val < 0.5 else 0)
        elif prop == 'Caco2_Wang':
            substituted.append(1 if val > -5.150 else 0)
        elif prop in ['CYP2C9_Substrate_CarbonMangels', 'CYP2D6_Substrate_CarbonMangels', 'CYP3A4_Substrate_CarbonMangels']:
            substituted.append(1 if 0.5 <= val <= 1 else 0)
    return substituted

# === Define ADMET-score calculation function ===
def calculate_admet_score(row):
    substituted = substitute_values(row)
    weighted_sum = sum(val * weights[i] for i, val in enumerate(substituted))
    return weighted_sum / sum(weights), substituted

# === Apply the scoring ===
results = df.apply(calculate_admet_score, axis=1)
df['ADMET_Score'] = results.apply(lambda x: x[0])
df_substituted = pd.DataFrame(results.apply(lambda x: x[1]).tolist(), columns=property_names)

# === Combine results and save ===
final_df = pd.concat([df[['SMILES', 'ADMET_Score']], df_substituted], axis=1)
final_df.to_csv(output_csv, index=False)

print(f"✅ ADMET-scores calculated and saved to: {output_csv}")