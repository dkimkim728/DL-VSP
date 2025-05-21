# Multi-Omics-Guided Deep Learning Pipeline for Inhibitor Discovery

This repository provides a comprehensive protocol for conducting **virtual drug screening** using a **multi-omics-integrated AI pipeline**. The framework combines **pharmacotranscriptomic regression**, **structure-based docking**, **molecular dynamics (MD) simulations**, and **ADMET filtering** to prioritize novel small-molecule inhibitors for experimental validation.

---

## 🧠 Key Features

- End-to-end pipeline integrating transcriptomic, structural, and pharmacokinetic data  
- Fine-tuned **MolFormer** model for molecular property prediction  
- High-throughput screening of chemical libraries (e.g., ZINC20, PubChem)  
- Robust post-processing with **docking**, **MD simulations**, and **ADMET scoring**  
- Easily adaptable for different target proteins and disease contexts

---

## 🏗️ Pipeline Overview

### 1. Pharmacotranscriptomic Screening
- **Input**: Gene expression profiles (e.g., LINCS L1000)  
- **Output**: Filtered gene–compound tables for target protein expression

### 2. Model Training
- **Stages**: Full dataset → Biological subset → Target-specific subset  
- **Model**: Hierarchically fine-tuned MolFormer  
- **Output**: Model checkpoints and prediction script

### 3. Structure-Based Virtual Screening
- **Input**: Target protein structure (PDB) and candidate ligand SMILES  
- **Tool**: AutoDock Vina  
- **Output**: Docking score (binding affinity)

### 4. Molecular Dynamics Simulation
- **Tool**: OpenMM  
- **Output**: Structural stability (RMSD)

### 5. ADMET Property Evaluation
- **Tool**: ADMET-AI  
- **Output**: Multi-property drug-likeness score

### 6. Prioritized Compound Ranking
- **Input**: Combined transcriptomic, docking, MD, and ADMET scores  
- **Output**: Ranked list of candidate inhibitors

---

## 🧬 Supported Data Sources

- **Transcriptomics**: LINCS CLUE platform, GEO, or custom datasets  
- **Chemicals**: ZINC20, PubChem, or in-house compound libraries  
- **Structures**: PDB (target protein), RDKit for ligand generation

---

## 🔧 Dependencies

- Python ≥ 3.8  
- [RDKit](https://www.rdkit.org)  
- [MolFormer](https://ibm.box.com/v/MolFormer-data)  
- [AutoDock Vina](https://vina.scripps.edu)  
- [Open Babel](http://openbabel.org)  
- [OpenMM](https://openmm.org)  
- [ADMET-AI](https://admet.ai.greenstonebio.com)  
- PyMOL / ChimeraX (for structural visualization)

> Create the environment using:  
> `conda env create -f env.yaml`

---

## 📂 Directory Structure

```plaintext
project_root/
│
├── 1_data_processing/
│   ├── parse_gctx_to_csv.py
│   ├── preprocess_expression_data.py
│   ├── calculate_logfc.py
│   └── outputs/
│
├── 2_model_training/
│   └── finetuning/
│       ├── finetune.py
│       ├── validate.py
│       ├── predict.py
│       ├── args.py
│       ├── utils.py
│       ├── tokenizer.py
│       ├── trainer.py
│       ├── run_finetune.sh
│       ├── run_validate_mu.sh
│       ├── run_predict_mu.sh
│       │
│       ├── rotate_attention/
│       │   ├── attention_layer.py
│       │   ├── rotary.py
│       │   ├── rotate_builder.py
│       │
│       ├── pretrained/
│       │   └── N-Step-Checkpoint_3_30000.ckpt
│       │
│       ├── source/
│           ├── sample_stage1.csv
│           ├── sample_stage2.csv
│           ├── sample_stage3.csv
│           ├── sample_val.csv
│           ├── sample_test.csv
│
├── 3_smiles_preparation/
│   ├── screen_valid_smiles.py
│   ├── generate_multiconformers.py
│   └── outputs/
│
├── 4_docking/
│   ├── run_docking_batch.py
│   └── outputs/
│
├── 5_md_simulation/
│   ├── convert_pdb_to_mol2.py
│   ├── mol2_to_frcmod.py
│   ├── assemble_complex_tleap.py
│   ├── run_md_simulation.py
│   ├── analyze_rmsd.py
│   ├── tleap.in
│   └── outputs/
│
├── 6_admet_scoring/
│   ├── calculate_admet_score.py
│   └── outputs/
│
├── data/
├── envs/
│   └── env.yaml
├── results/
└── .gitignore
