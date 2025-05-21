"""
run_md_simulation.py

Runs molecular dynamics using OpenMM from AMBER input files (.prmtop and .inpcrd).
Outputs: trajectory (.dcd), scalar energies (.csv), minimized structure (.pdb).
"""

import os
from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout

# === Input Configuration ===
complex_id = "ligand_complex"  # <-- Update to match your folder name

root_dir = os.path.join("5_md_simulation", "inputs", complex_id)
output_dir = os.path.join("5_md_simulation", "outputs", complex_id)

inpcrd_file = os.path.join(root_dir, "ligand_wat.inpcrd")
prmtop_file = os.path.join(root_dir, "ligand_wat.prmtop")
init_pdb_path = os.path.join(output_dir, "init.pdb")
traj_path = os.path.join(output_dir, "traj.dcd")
scalars_path = os.path.join(output_dir, "scalars.csv")

# === Ensure output directory exists ===
os.makedirs(output_dir, exist_ok=True)

# === Load input system ===
inpcrd = AmberInpcrdFile(inpcrd_file)
prmtop = AmberPrmtopFile(prmtop_file, periodicBoxVectors=inpcrd.boxVectors)
system = prmtop.createSystem(
    nonbondedMethod=PME,
    nonbondedCutoff=1*nanometer,
    constraints=HBonds
)
system.addForce(MonteCarloBarostat(1*bar, 300*kelvin))

# === Setup Integrator and Simulation ===
integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.004*picoseconds)
simulation = Simulation(prmtop.topology, system, integrator)
simulation.context.setPositions(inpcrd.positions)

# === Minimize Energy and Save Initial Structure ===
print(f"🔬 Minimizing energy for {complex_id}...")
simulation.minimizeEnergy()
positions = simulation.context.getState(getPositions=True).getPositions()
with open(init_pdb_path, 'w') as f:
    PDBFile.writeFile(simulation.topology, positions, f)
print(f"✅ Minimized structure saved to: {init_pdb_path}")

# === Attach Reporters ===
simulation.reporters.append(DCDReporter(traj_path, 10))
simulation.reporters.append(StateDataReporter(stdout, 100, step=True, temperature=True, elapsedTime=True))
simulation.reporters.append(StateDataReporter(
    scalars_path, 10,
    step=True, time=True,
    potentialEnergy=True, totalEnergy=True,
    temperature=True
))

# === Run MD Simulation ===
print(f"🚀 Running simulation for 1,000,000 steps...")
simulation.step(1_000_000)
print(f"✅ Simulation completed. Results in: {output_dir}")