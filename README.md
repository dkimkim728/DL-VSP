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
- **Model**: Fine-tuned MolFormer  
- **Output**: Transcriptomic suppression score (e.g., LFC)

### 2. Model Training
- **Stages**: Full dataset → Biological subset → Target-specific subset  
- **Tool**: Hierarchically fine-tuned MolFormer  
- **Output**: Model checkpoints and prediction script

### 3. Library Input
- **Input**: Canonical SMILES for candidate ligands  
- **Tools**: RDKit, Open Babel  
- **Output**: Validated, conformer-enriched ligand files

### 4. Structure-Based Virtual Screening
- **Input**: Target protein (PDBQT) and ligand library (PDBQT)  
- **Tool**: AutoDock Vina  
- **Output**: Docking scores (binding affinity)

### 5. Molecular Dynamics Simulation
- **Input**: Docked complexes  
- **Tool**: OpenMM  
- **Output**: Trajectory (DCD), minimized PDB, energy profiles

### 6. ADMET Property Evaluation
- **Input**: Ligand library  
- **Tool**: ADMET-AI  
- **Output**: Drug-likeness score across 41 features

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
> `conda env create -f envs/env.yaml`

---

## 📂 Directory Structure

```plaintext
Protocol/
│
├── 1_data_processing/
│   ├── calculate_logfc.py
│   ├── parse_gctx_to_csv.py
│   └── preprocess_expression_data.py
│
├── 2_model_training/
│   ├── finetune.py
│   ├── validate.py
│   ├── predict.py
│   ├── run_finetune.sh
│   ├── run_validate.sh
│   ├── run_predict.sh
│   ├── args.py
│   ├── tokenizer.py
│   ├── trainer.py
│   ├── utils.py
│   └── finetuning/
│       ├── pretrained/
│       │   └── N-Step-Checkpoint_3_30000.ckpt
│       ├── rotate_attention/
│       │   ├── attention_layer.py
│       │   ├── rotary.py
│       │   └── rotate_builder.py
│       └── source/
│           ├── sample_stage1.csv
│           ├── sample_stage2.csv
│           ├── sample_stage3.csv
│           ├── sample_val.csv
│           └── sample_test.csv
│
├── 3_library_input/
│   ├── generate_multiconformers.py
│   ├── predict.py
│   └── screen_valid_smiles.py
│
├── 4_docking/
│   ├── run_docking_batch.py
│   ├── input/sample/
│   │   ├── conf.txt
│   │   ├── sample_ligand.pdb
│   │   └── sample_protein.pdb
│   └── output/sample/
│       ├── sample_variant0_conf0.pdbqt
│       ├── sample_variant0_conf1.pdbqt
│       ├── ...
│       └── sample_variant0_conf9.pdbqt
│
├── 5_md_simulation/
│   ├── analyze_rmsd.py
│   ├── convert_pdb_to_mol2.py
│   ├── mol2_to_frcmod.py
│   ├── run_md_simulation.py
│   ├── tleap.in
│   ├── input/sample/
│   │   ├── frcmod/ligand.frcmod
│   │   ├── gaff_mol2/ligand.mol2
│   │   ├── ligand_complex/
│   │   │   ├── ligand_vac.pdb
│   │   │   ├── ligand_wat.inpcrd
│   │   │   ├── ligand_wat.pdb
│   │   │   └── ligand_wat.prmtop
│   │   └── pdb/
│   │       ├── sample_ligand.pdb
│   │       └── sample_protein.pdb
│   └── output/sample/
│       ├── init.pdb
│       ├── leap.log
│       ├── RMSD_ligand.png
│       └── scalars.csv
│
├── 6_admet_scoring/
│   ├── calculate_admet_score.py
│   ├── input/sample/
│   │   └── smiles.csv
│   └── output/sample/
│       ├── admet_score.csv
│       └── result.csv
│
├── data/
├── envs/
│   └── env.yaml
├── results/
├── LICENSE
└── .gitignore
