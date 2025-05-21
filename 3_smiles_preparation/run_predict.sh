#!/bin/bash

# Inference script using trained MolFormer model for large-scale prediction on SMILES

python predict.py \
  --seed_path "./checkpoints/stage_3_full_model.pt" \
  --input_file "./data/library.csv" \
  --output_file "./results/library_predictions.csv" \
  --measure_name "FC" \
  --batch_size 512 \
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
  --checkpoint_root "./checkpoints" \
  --data_root "./data"