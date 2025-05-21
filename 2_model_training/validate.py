"""
validate.py

Evaluate a fine-tuned MolFormer model on validation or test datasets.
Outputs MAE, RMSE, average loss, and residuals to CSV.
"""

import os
import torch
import numpy as np
import pandas as pd
from trainer import parse_args, PropertyPredictionDataModule
from sklearn.metrics import mean_absolute_error, mean_squared_error

def evaluate_model(model, data_loader, device):
    model.eval()
    all_targets = []
    all_predictions = []
    total_loss = 0.0
    loss_fn = model.loss_function
    num_samples = 0

    with torch.no_grad():
        for batch in data_loader:
            inputs, attention_mask, targets = batch
            inputs = inputs.to(device)
            attention_mask = attention_mask.to(device)
            targets = targets.to(device)

            outputs = model(inputs)
            preds = outputs.mean(dim=1)

            loss = loss_fn(preds, targets)
            total_loss += loss.item() * inputs.size(0)

            all_predictions.extend(preds.cpu().numpy())
            all_targets.extend(targets.cpu().numpy())
            num_samples += inputs.size(0)

    avg_loss = total_loss / num_samples
    mae = mean_absolute_error(all_targets, all_predictions)
    rmse = np.sqrt(mean_squared_error(all_targets, all_predictions))

    return avg_loss, mae, rmse, all_targets, all_predictions

def main():
    args = parse_args()
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"📟 Running on: {device}")

    # Load data
    datamodule = PropertyPredictionDataModule(args)
    datamodule.prepare_data()
    val_loader = datamodule.val_dataloader()[0]

    # Load model
    model_path = os.path.join(args.checkpoint_root, "stage_3_full_model.pt")
    if not os.path.exists(model_path):
        raise FileNotFoundError(f"🚨 Model checkpoint not found at: {model_path}")
    
    print(f"✅ Loading model from: {model_path}")
    model = torch.load(model_path, map_location=device)
    model.to(device)

    # Run evaluation
    print("🔍 Evaluating...")
    avg_loss, mae, rmse, targets, predictions = evaluate_model(model, val_loader, device)

    print(f"\n📊 Evaluation Results:")
    print(f"  • Average Loss: {avg_loss:.4f}")
    print(f"  • MAE: {mae:.4f}")
    print(f"  • RMSE: {rmse:.4f}")

    # Save residuals to CSV
    residuals = np.array(predictions) - np.array(targets)
    df_out = pd.DataFrame({
        "Actual": targets,
        "Predicted": predictions,
        "Residual": residuals
    })

    output_csv = os.path.join(args.checkpoint_root, "residuals_validation.csv")
    df_out.to_csv(output_csv, index=False)
    print(f"\n📝 Residuals saved to: {output_csv}")

if __name__ == '__main__':
    main()