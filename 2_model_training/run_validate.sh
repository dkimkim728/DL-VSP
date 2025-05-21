#!/bin/bash

# Validation script for the final fine-tuned MolFormer model
# Outputs MAE, RMSE, and residual summary from validation set

python validate.py \
  --valid_data_path "source/sample_val.csv" \
  --measure_name "FC" \
  --dataset_name "multi_stage_validation" \
  --batch_size 64 \
  --n_embd 512 \
  --n_head 8 \
  --n_layer 6 \
  --num_feats 128 \
  --dropout 0.2 \
  --d_dropout 0.1 \
  --lr_start 1e-4 \
  --lr_multiplier 1.0 \
  --max_epochs 5 \
  --num_workers 4 \
  --checkpoint_root "checkpoints/" \
  --data_root "source/"