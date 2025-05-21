"""
analyze_rmsd.py

Compares RMSD trajectories between a ligand-protein complex and an apo-protein system.
Requires pre-generated MD trajectories (.dcd) and topologies (.prmtop).
"""

import os
import mdtraj as md
import matplotlib.pyplot as plt

# === Configurable Input Paths ===
traj_dir = "5_md_simulation/inputs/"
topo_dir = "5_md_simulation/outputs/"
output_plot = "5_md_simulation/outputs/rmsd_comparison.png"

# Complex system
complex_name = "complex"  # e.g., com289
complex_traj = os.path.join(traj_dir, complex_name, "traj.dcd")
complex_top = os.path.join(topo_dir, complex_name, "ligand_wat.prmtop")

# Protein-only system
protein_name = "apo"  # e.g., onlyP
protein_traj = os.path.join(traj_dir, protein_name, "traj.dcd")
protein_top = os.path.join(topo_dir, protein_name, "protein_wat.prmtop")

# === Load and process complex trajectory ===
top_complex = md.load_prmtop(complex_top)
ref_complex = md.load(complex_traj, top=top_complex, frame=0)
sel_complex = ref_complex.topology.select("protein")
ref_sel_complex = ref_complex.atom_slice(sel_complex)

rmsds_complex = []
for chunk in md.iterload(complex_traj, top=top_complex, chunk=100):
    sel_chunk = chunk.atom_slice(sel_complex)
    rmsd_chunk = md.rmsd(sel_chunk, ref_sel_complex)
    rmsds_complex.extend(rmsd_chunk)

# === Load and process apo-protein trajectory ===
top_protein = md.load_prmtop(protein_top)
ref_protein = md.load(protein_traj, top=top_protein, frame=0)
sel_protein = ref_protein.topology.select("protein")
ref_sel_protein = ref_protein.atom_slice(sel_protein)

rmsds_protein = []
for chunk in md.iterload(protein_traj, top=top_protein, chunk=100):
    sel_chunk = chunk.atom_slice(sel_protein)
    rmsd_chunk = md.rmsd(sel_chunk, ref_sel_protein)
    rmsds_protein.extend(rmsd_chunk)

# === Plotting ===
plt.figure(figsize=(10, 7))
plt.plot(rmsds_protein, label="Apo-protein", color="black", linewidth=2)
plt.plot(rmsds_complex, label="Complex", color="blue", linewidth=2)
plt.xlabel("Time (ps)", fontsize=13)
plt.ylabel("RMSD (nm)", fontsize=13)
plt.title("RMSD Comparison", fontsize=14)
plt.legend(fontsize=12)
plt.grid(True)
plt.tight_layout()
plt.savefig(output_plot, dpi=300)
plt.show()

print(f"✅ RMSD comparison plot saved to: {output_plot}")