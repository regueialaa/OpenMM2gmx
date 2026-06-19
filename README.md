# OpenMM2gmx

<img width="1024" height="559" alt="image" align="center" src="https://github.com/user-attachments/assets/0f7ccc06-e66b-4ca1-91b5-4342222cc7e3" />

A workflow for converting and centering OpenMM trajectories using GROMACS.
Works for any system and has been tested on membrane proteins with a triclinic box.

## Overview

This repository provides a workflow to:

- Generate `.gro` , `.top ` and `.tpr` files from OpenMM outputs
- Convert OpenMM `.dcd` trajectories into GROMACS-compatible `.xtc` , handling esapcially difference between triclinic box definiton between openmm and gromcas
- treating the the trajectory with no jump and center options around a selection of residues (eg. Center of mass of the protein)

## Author

Alaa REGUEI, PhD Student - Université Paris Cité, BFA.

## Installation

```bash
pip install -r requirements.txt
```
also you will need to install parmed and openmm ,we recommend using conda ;

Make sure you have GROMACS installed and available in your `PATH` (e.g. `gmx trjconv`, `gmx make_ndx`).
To install it  , please refer to official guide : https://manual.gromacs.org/current/install-guide/index.html

## Usage

Run the script from the command line with the following arguments:

```bash
python OpenMM2gmx.py \
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

## Command-Line Arguments

--top (required)
Input topology file (e.g. parm7)

--xml (required)
OpenMM coordinates/system XML file from the production run

--traj (required)
One or more OpenMM trajectory files (e.g. traj.dcd or traj_part1.dcd traj_part2.dcd ...). You can pass multiple files.

--mdp (required)
Input GROMACS .mdp file

--center_res (required)
Residue(s) used to center the simulation. Use COM to center on the protein center of mass. Use residue ranges like 10-50, single residues like 10, or multiple residues/ranges like 10 20 30 or 10-20 30-40. Examples: --center_res COM, --center_res 10-50, --center_res 10 20 30

--sim_name (required)
Name of the simulation/replica (used in output filenames and logs)

--stride (optional)
Stride for trajectory analysis (default: 1). Not recommended if the trajectory is partitioned.

--save_mode (optional)
Output mode (default: 0). 0 = full corrected system trajectory (Protein + non-polymer), 1 = protein-only corrected trajectory, 2 = both trajectories

--output_dir (optional)
Output directory name (default: MD_output)
