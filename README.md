# OpenMM2gmx

<!-- <img width="1024" height="559" alt="image" align="center" src="https://github.com/user-attachments/assets/0f7ccc06-e66b-4ca1-91b5-4342222cc7e3" /> -->
<img width="1024" height="436" alt="image" align="center" src="https://github.com/user-attachments/assets/bc48fe8c-4316-4f39-a4d5-4593d9c21176" />


A workflow for converting and centering OpenMM trajectories using GROMACS.
Works for any system and has been tested on membrane proteins with a triclinic box.

## Overview

This repository provides a workflow to:

- Generate `.gro` , `.top ` and `.tpr` files from OpenMM outputs
- Convert OpenMM `.dcd` trajectories into GROMACS-compatible `.xtc` , handling esapcially difference between triclinic box definiton between openmm and gromcas
- treating the the trajectory with no jump and center options around a selection of residues (eg. Center of mass of the protein)

## Authors

Alaa REGUEI, PhD Student - Université Paris Cité, BFA.\
Samuel Murail, Associate Professor - Université Paris Cité, BFA.

## Installation

You will need first to install **parmed** and **openmm** ,we recommend using conda ;

Make sure you have **GROMACS** installed and available in your `PATH` (e.g. `gmx trjconv`, `gmx make_ndx`).
To install it  , please refer to official guide : https://manual.gromacs.org/current/install-guide/index.html

```bash
git clone https://github.com/regueialaa/OpenMM2gmx.git
cd OpenMM2gmx
pip install .
OpenMM2gmx --help
```

## Usage

OpenMM2gmx can be used either from the command line or directly in Python (e.g., in a Jupyter Notebook).

### Option 1: Command line

Run the script with the following arguments:

```bash
OpenMM2gmx \
  --top topfile.parm7 \
  --xml system.xml \
  --traj traj1.dcd traj2.dcd \
  --mdp simulation.mdp \
  --center_res COM \
  --sim_name my_simulation \
  --stride 1 \
  --save_mode 0 \
  --output_dir MD_output
```

### Option 2: Python / Jupyter Notebook

OpenMM2gmx can also be imported and used directly in Python scripts or Jupyter Notebooks.

An example notebook demonstrating the package on a Class B GPCR system is provided:

```text
GPCR_case.ipynb
```

This tutorial shows how to import the package and run the workflow directly from a Python environment.

## Arguments

--top (required)
Input topology file : eg. parm7 (recommanded for membrane systems) , pdb , cif or mmcif

--xml (required)
OpenMM coordinates/system XML file from the production run

--traj (required)
One or more OpenMM trajectory files (e.g. traj.dcd or traj_part1.dcd traj_part2.dcd ...). You can pass multiple files.

--mdp (required)
Input GROMACS .mdp file

--center_res (required)
Residue used to define the simulation center: "COM" centers on the residue closest to the center of mass of the system, or you can specify your own residue in the form "resid100". Example: --center_res "COM" (default) or --center_res "resid:A:100"'

--sim_name (required)
Name of the simulation/replica (used in output filenames and logs)

--stride (optional)
Stride for trajectory analysis (default: 1). Not recommended if the trajectory is partitioned.

--save_mode (optional)
Output mode (default: 0). 0 = full corrected system trajectory (Protein + non-polymer), 1 = protein-only corrected trajectory, 2 = both trajectories

--output_dir (optional)
Output directory name (default: MD_output)

## References

If you use OpenMM2gmx in your work, please cite the following software packages:

- **OpenMM**: Eastman *et al.* (2024). *OpenMM 8: Molecular Dynamics Simulation with Machine Learning Potentials*. The Journal of Physical Chemistry B, 128 (1), 109–116.

- **GROMACS**: Abraham *et al.* (2015). *GROMACS: High performance molecular simulations through multi-level parallelism from laptops to supercomputers*. SoftwareX, 1-2, 19–25.
