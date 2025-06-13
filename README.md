# Deep Learning-Enabled Virtual Drug Screening Protocol (DL-VSP)

This repository provides a reproducible framework for conducting **virtual drug screening** using a **deep learning–enabled pipeline for protein-targeted therapeutics**. The DL-VSP protocol integrates **chemical language modeling**, **structure-based docking**, **molecular dynamics (MD) simulations**, and **ADMET filtering** to prioritize candidate small-molecule inhibitors from ultra-large compound libraries.

The pipeline is generalizable across disease areas and compatible with canonical SMILES-based input formats. A fine-tuned MolFormer model is used here as an example predictive engine for regression tasks such as transcriptomic modulation.

---

## 🧠 Key Features

- End-to-end screening pipeline from chemical structure to candidate selection  
- Integration of **MolFormer**, **AutoDock Vina**, **OpenMM**, and **ADMET-AI**  
- Hierarchical model fine-tuning on multi-scale datasets (general to disease-specific)  
- Compatible with ultra-large libraries such as ZINC20  
- Modular design supporting diverse bioactivity prediction tasks

---

## 🏗️ Pipeline Overview

1. **Data Curation**  
   - Input: Transcriptomic or bioactivity datasets (e.g., LINCS)  
   - Output: SMILES-labeled CSVs with standardized activity values (e.g., LFC)  

2. **Model Fine-Tuning**  
   - Tool: MolFormer pretrained on 100M compounds (ZINC + PubChem)  
   - Phases: General → Disease-specific → Target-specific  
   - Output: Fine-tuned model checkpoint for regression/classification tasks  

3. **Compound Library Preparation**  
   - Source: ZINC20 (https://files.docking.org/zinc20-ML/)  
   - Tools: `screen_valid_smiles.py`, `generate_multiconformers.py`  
   - Output: Validated and expanded SMILES inputs (PDB, PDBQT)

4. **Structure-Based Docking**  
   - Tool: AutoDock Vina (selected for scalability)  
   - Output: Docking scores per conformer (BA in kcal/mol)

5. **Molecular Dynamics Simulation**  
   - Tool: OpenMM + AMBER-style prep (via tleap)  
   - Output: Time-resolved RMSD and trajectory stability metrics

6. **ADMET Filtering**  
   - Tool: ADMET-AI (GNN-based)  
   - Output: 41 pharmacokinetic and toxicity predictions, including composite ADMET scores

---

## 📊 Representative Results

The DL-VSP pipeline was benchmarked on a multi-stage fine-tuning task targeting **PIN1 inhibition in TNBC**:

- **Figure 3**: Multi-phase MAE loss curves and residual distributions (Pancancer → BC → TNBC stages) demonstrate progressive refinement.
- **Figure 4**: Compound reduction funnel shows stepwise filtering across model prediction (LFC < –2.0), docking (BA < –6.0), MD (RMSD < 0.6 nm), and ADMET (score > 0.5), yielding 7 prioritized candidates from an initial 10M ZINC20 entries.

---

## 🧬 Supported Input Sources

- **Expression datasets**: LINCS, GEO, DepMap, ENCODE  
- **Compounds**: ZINC20, PubChem, in-house libraries  
- **Structure files**: PDB, MOL2, or SMILES

---

## 🔧 Environment Setup

Each module operates within its own Conda environment:

```bash
conda env create -f envs/MolFormer_env.yaml     # Model fine-tuning
conda env create -f envs/Vina_env.yaml          # Docking
conda env create -f envs/Ambertools_env.yaml    # MD simulation
conda env create -f envs/ADMETai_env.yaml       # ADMET prediction
```

Activate with `conda activate <env_name>`.

---

## 📂 Repository Structure

\`\`\`plaintext
Protocol/
├── 1_data_processing/              # Expression data parsing and labeling
├── 2_model_training/              # MolFormer fine-tuning (multi-phase)
├── 3_library_input/               # SMILES standardization, conformer generation
├── 4_protein_ligand_docking/      # Docking scripts and config
├── 5_md_simulation/               # Force field setup, simulation, RMSD analysis
├── 6_admet_scoring/               # ADMET-AI prediction and scoring
├── data/                          # Raw data and preprocessed samples
├── envs/                          # Conda environment YAMLs
├── results/                       # Output predictions and filtered candidates
└── LICENSE, .gitignore
\`\`\`

---

## 🔬 Reference Model: MolFormer

The pretrained MolFormer model was developed by IBM and trained on ~100M canonical SMILES from ZINC and PubChem. It uses linear attention with rotary embeddings and supports downstream fine-tuning for both regression and classification tasks.

Checkpoints available at: https://ibm.box.com/v/MolFormer-data  
Original study: Ross *et al.* (2022), *Nature Machine Intelligence*  
DOI: [10.1038/s42256-022-00580-7](https://doi.org/10.1038/s42256-022-00580-7)

---

## 📄 Citation

If you use this repository, please cite the associated MoLFormer paper and acknowledge this pipeline in your work. See [CITATION.cff](./CITATION.cff) for formatted references.
