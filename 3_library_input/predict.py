"""
predict.py

Predict regression values for a large SMILES library using the final fine-tuned model.
Outputs a CSV file with predicted values for each SMILES entry.
"""

import torch
import pandas as pd
from trainer import parse_args, tokenize_smiles

def predict_smiles_from_csv(model, tokenizer, input_csv, smiles_column, output_csv, device, batch_size=512):
    df = pd.read_csv(input_csv)
    smiles_data = df[smiles_column].dropna().astype(str).tolist()
    
    if len(smiles_data) == 0:
        raise ValueError(f"No valid SMILES strings found in column '{smiles_column}'.")

    model.eval()
    predictions = []

    with torch.no_grad():
        for i in range(0, len(smiles_data), batch_size):
            batch_smiles = smiles_data[i:i + batch_size]
            input_ids, attention_mask = tokenize_smiles(batch_smiles, tokenizer, device)
            torch.cuda.empty_cache()

            batch_preds = model(input_ids).mean(dim=1).cpu().numpy().flatten()
            predictions.extend(batch_preds)

            print(f"✅ Processed {i + len(batch_smiles)}/{len(smiles_data)} molecules")

    output_df = pd.DataFrame({
        smiles_column: smiles_data,
        "Predicted_LFC": predictions
    })
    output_df.to_csv(output_csv, index=False)
    print(f"✅ Saved predictions to: {output_csv}")

def main():
    args = parse_args()
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    # Load model
    model = torch.load(args.seed_path, map_location=device)
    model.eval()

    # Predict from input file
    predict_smiles_from_csv(
        model=model,
        tokenizer=model.tokenizer,
        input_csv=args.input_file,
        smiles_column="canonical_smiles",
        output_csv=args.output_file,
        device=device
    )

if __name__ == '__main__':
    main()