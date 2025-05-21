import os
import torch
import numpy as np
import pandas as pd
import regex as re
import logging
from torch import nn
from torch.utils.data import DataLoader, Dataset
from sklearn.metrics import mean_absolute_error, mean_squared_error
from transformers import BertTokenizer
import pytorch_lightning as pl
from fast_transformers.masking import LengthMask as LM
from rotate_attention.rotate_builder import RotateEncoderBuilder
from fast_transformers.feature_maps import GeneralizedRandomFeatures
from functools import partial
import argparse
from rdkit import Chem

# ==================== Tokenizer ====================

PATTERN = "(\[[^\]]+]|Br?|Cl?|N|O|S|P|F|I|b|c|n|o|s|p|\(|\)|\.|=|#|-|\+|\\\\|\/|:|~|@|\?|>|\*|\$|\%[0-9]{2}|[0-9])"

class MolTranBertTokenizer(BertTokenizer):
    def __init__(self, vocab_file="bert_vocab.txt", **kwargs):
        super().__init__(
            vocab_file,
            unk_token="<pad>",
            sep_token="<eos>",
            pad_token="<pad>",
            cls_token="<bos>",
            mask_token="<mask>",
            **kwargs
        )
        self.regex_tokenizer = re.compile(PATTERN)
        self.wordpiece_tokenizer = None
        self.basic_tokenizer = None

    def _tokenize(self, text):
        return self.regex_tokenizer.findall(text)

    def convert_tokens_to_string(self, tokens):
        return "".join(tokens).strip()

# ==================== Utility ====================

def normalize_smiles(smiles, canonical=True, isomeric=False):
    try:
        mol = Chem.MolFromSmiles(smiles)
        return Chem.MolToSmiles(mol, isomericSmiles=isomeric) if mol else None
    except:
        return None

# ==================== Dataset ====================

class PropertyPredictionDataset(Dataset):
    def __init__(self, df, measure_name, tokenizer, aug=True):
        df = df.dropna(subset=['canonical_smiles', measure_name])
        self.original_smiles = df["canonical_smiles"].tolist()
        self.original_canonical_map = {
            smi: normalize_smiles(smi, canonical=True, isomeric=False) for smi in self.original_smiles
        }
        self.measure_map = dict(zip(self.original_smiles, df[measure_name].tolist()))
        self.tokenizer = tokenizer
        self.aug = aug

    def __getitem__(self, index):
        original_smiles = self.original_smiles[index]
        canonical_smiles = self.original_canonical_map[original_smiles]
        return canonical_smiles, self.measure_map[original_smiles]

    def __len__(self):
        return len(self.original_smiles)

class PropertyPredictionDataModule(pl.LightningDataModule):
    def __init__(self, hparams):
        super().__init__()
        self.hparams = hparams
        self.tokenizer = MolTranBertTokenizer("bert_vocab.txt")
        self.measure_name = hparams.measure_name

    def prepare_data(self):
        self.train_ds = PropertyPredictionDataset(
            pd.read_csv(self.hparams.train_data_path),
            self.measure_name,
            self.tokenizer,
            aug=self.hparams.aug
        )
        self.val_ds = PropertyPredictionDataset(
            pd.read_csv(self.hparams.valid_data_path),
            self.measure_name,
            self.tokenizer,
            aug=False
        )
        self.test_ds = PropertyPredictionDataset(
            pd.read_csv(self.hparams.test_data_path),
            self.measure_name,
            self.tokenizer,
            aug=False
        )

    def collate(self, batch):
        valid_batch = [(s, t) for s, t in batch if isinstance(s, str)]
        if not valid_batch:
            return None
        smiles, targets = zip(*valid_batch)
        tokens = self.tokenizer.batch_encode_plus(smiles, padding=True, add_special_tokens=True)
        return (
            torch.tensor(tokens['input_ids'], dtype=torch.long),
            torch.tensor(tokens['attention_mask'], dtype=torch.long),
            torch.tensor(targets, dtype=torch.float)
        )

    def train_dataloader(self):
        return DataLoader(self.train_ds, batch_size=self.hparams.batch_size, shuffle=True, num_workers=self.hparams.num_workers, collate_fn=self.collate)

    def val_dataloader(self):
        return DataLoader(self.val_ds, batch_size=self.hparams.batch_size, shuffle=False, num_workers=self.hparams.num_workers, collate_fn=self.collate)

# ==================== Model ====================

class LightningModule(pl.LightningModule):
    def __init__(self, config, tokenizer):
        super().__init__()
        self.save_hyperparameters(config)
        n_vocab, d_emb = len(tokenizer.vocab), config.n_embd

        builder = RotateEncoderBuilder.from_kwargs(
            n_layers=config.n_layer,
            n_heads=config.n_head,
            query_dimensions=config.n_embd // config.n_head,
            value_dimensions=config.n_embd // config.n_head,
            feed_forward_dimensions=config.n_embd,
            attention_type='linear',
            feature_map=partial(GeneralizedRandomFeatures, n_dims=config.num_feats),
            activation='gelu'
        )

        self.tok_emb = nn.Embedding(n_vocab, config.n_embd)
        self.drop = nn.Dropout(config.d_dropout)
        self.blocks = builder.get()
        self.net = nn.Sequential(
            nn.Linear(config.n_embd, config.n_embd),
            nn.GELU(),
            nn.Dropout(config.dropout),
            nn.Linear(config.n_embd, 1)
        )
        self.loss_fn = nn.L1Loss()

    def forward(self, inputs):
        x = self.tok_emb(inputs)
        x = self.drop(x)
        x = self.blocks(x)
        x = self.net(x).squeeze()
        return x

    def training_step(self, batch, batch_idx):
        x, mask, y = batch
        x = self(x)
        loss = self.loss_fn(x, y)
        return {"loss": loss}

    def validation_step(self, batch, batch_idx):
        x, mask, y = batch
        preds = self(x)
        loss = self.loss_fn(preds, y)
        return {"val_loss": loss}

    def configure_optimizers(self):
        return torch.optim.AdamW(self.parameters(), lr=self.hparams.lr_start)

# ==================== Argparser ====================

def parse_args():
    parser = argparse.ArgumentParser(description="MolFormer Fine-Tuning")
    parser.add_argument("--train_data_path", type=str, required=True)
    parser.add_argument("--train_data_path_stage2", type=str)
    parser.add_argument("--train_data_path_stage3", type=str)
    parser.add_argument("--valid_data_path", type=str, required=True)
    parser.add_argument("--test_data_path", type=str, required=True)
    parser.add_argument("--checkpoint_root", type=str, required=True)
    parser.add_argument("--dataset_name", type=str, default="molformer")
    parser.add_argument("--measure_name", type=str, required=True)
    parser.add_argument("--batch_size", type=int, default=32)
    parser.add_argument("--n_embd", type=int, default=768)
    parser.add_argument("--n_head", type=int, default=12)
    parser.add_argument("--n_layer", type=int, default=6)
    parser.add_argument("--num_feats", type=int, default=32)
    parser.add_argument("--d_dropout", type=float, default=0.1)
    parser.add_argument("--dropout", type=float, default=0.1)
    parser.add_argument("--lr_start", type=float, default=3e-5)
    parser.add_argument("--max_epochs", type=int, default=5)
    parser.add_argument("--num_workers", type=int, default=8)
    parser.add_argument("--device", type=str, default="cuda")
    parser.add_argument("--aug", type=bool, default=False)
    return parser.parse_args()