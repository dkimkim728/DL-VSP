#!/bin/bash

# Multi-Stage Fine-Tuning Script for MolFormer (Using Sampled Datasets)
# Stage 1 → General dataset (sample_stage1.csv)
# Stage 2 → Subset dataset (sample_stage2.csv)
# Stage 3 → Target-specific dataset (sample_stage3.csv)

nice -n 19 python finetune.py \
    --device cuda \
    --batch_size 32 \
    --n_head 12 \
    --n_layer 6 \
    --n_embd 768 \
    --d_dropout 0.1 \
    --dropout 0.1 \
    --lr_start 3e-5 \
    --num_workers 8 \
    --max_epochs 5 \
    --num_feats 32 \
    --checkpoint_every 100 \
    --seed_path "pretrained/N-Step-Checkpoint_3_30000.ckpt" \
    --dataset_name "cmap_data" \
    --measure_name "FC" \
    --checkpoints_folder "checkpoints/" \
    --checkpoint_root "checkpoints/" \
    --data_path "source/" \
    --train_data_path "source/sample_stage1.csv" \
    --train_data_path_stage2 "source/sample_stage2.csv" \
    --train_data_path_stage3 "source/sample_stage3.csv" \
    --valid_data_path "source/sample_val.csv" \
    --test_data_path "source/sample_test.csv" \
    --data_root "source/" \
    --train_dataset_length 10000 \
    --seed 42 \
    --aug True