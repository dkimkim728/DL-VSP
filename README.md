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

### 2. Structure-Based Virtual Screening
- **Input**: Target protein structure (PDB) and candidate ligand SMILES
- **Tool**: AutoDock Vina
- **Output**: Docking score (binding affinity)

### 3. Molecular Dynamics Simulation
- **Tool**: OpenMM
- **Output**: Structural stability (RMSD)

### 4. ADMET Property Evaluation
- **Tool**: ADMET-AI
- **Output**: Multi-property drug-likeness score

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
├── data/                   # Preprocessed gene expression, SMILES, etc.
├── model/                  # Fine-tuned MolFormer checkpoints
├── docking/                # Ligand structures and docking configs
├── md/                     # Molecular dynamics simulation inputs/outputs
├── admet/                  # ADMET scoring results
├── scripts/                # Pipeline automation and analysis scripts
└── results/                # Final prioritized compounds
