"""
finetune.py

Implements three-stage hierarchical fine-tuning of MolFormer using SMILES and regression targets.
Designed to be executed via run_finetune.sh.
"""

import os
import torch
from trainer import LightningModule, parse_args, PropertyPredictionDataModule
from pytorch_lightning import Trainer
from pytorch_lightning.callbacks import ModelCheckpoint
from pytorch_lightning.loggers import TensorBoardLogger

def run_stage(stage_num, args, train_path, checkpoint_load=None, checkpoint_save=None):
    print(f"\n🚀 Starting Stage {stage_num} Training on {train_path}")
    args.train_data_path = train_path

    datamodule = PropertyPredictionDataModule(args)
    datamodule.prepare_data()

    if checkpoint_load and os.path.exists(checkpoint_load):
        print(f"📥 Loading checkpoint from {checkpoint_load}")
        model = torch.load(checkpoint_load, map_location=args.device)
    else:
        model = LightningModule(args, datamodule.tokenizer)

    checkpoint_callback = ModelCheckpoint(
        dirpath=os.path.dirname(checkpoint_save),
        filename=os.path.basename(checkpoint_save).replace(".pt", ""),
        save_top_k=1,
        monitor="val_loss",
        mode="min"
    )

    trainer = Trainer(
        max_epochs=args.max_epochs,
        accelerator="gpu" if torch.cuda.is_available() else "cpu",
        devices=1,
        callbacks=[checkpoint_callback],
        logger=TensorBoardLogger(save_dir=args.checkpoint_root, name=f"stage{stage_num}_logs")
    )

    trainer.fit(model, datamodule=datamodule)

    print(f"💾 Saving model to {checkpoint_save}")
    torch.save(model, checkpoint_save)

def main():
    args = parse_args()
    args.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    os.makedirs(args.checkpoint_root, exist_ok=True)

    stage1_ckpt = os.path.join(args.checkpoint_root, "stage_1_model.pt")
    stage2_ckpt = os.path.join(args.checkpoint_root, "stage_2_model.pt")
    stage3_ckpt = os.path.join(args.checkpoint_root, "stage_3_full_model.pt")

    run_stage(1, args, args.train_data_path, checkpoint_save=stage1_ckpt)
    run_stage(2, args, args.train_data_path_stage2, checkpoint_load=stage1_ckpt, checkpoint_save=stage2_ckpt)
    run_stage(3, args, args.train_data_path_stage3, checkpoint_load=stage2_ckpt, checkpoint_save=stage3_ckpt)

if __name__ == '__main__':
    main()