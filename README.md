# OpenMM2gmx

<img width="1024" height="559" alt="image" align="center" src="https://github.com/user-attachments/assets/0f7ccc06-e66b-4ca1-91b5-4342222cc7e3" />

A workflow for converting and centering OpenMM trajectories with GROMACS.
Adapted for membrane protein , run with triclinic box

---

## Overview

This repository provides a workflow to:

- Generate `.gro` , `.top ` and `.tpr` files from OpenMM outputs
- Convert OpenMM `.dcd` trajectories into GROMACS-compatible `.xtc` , handling esapcially difference between triclinic box definiton between openmm and gromcas
- treating the the trajectory with no jump and center options around a residue close to system center of mass or any other residues
- Visualize the final trajectory in VMD

---

## 1. Generate `.gro` , `.top ` and `.tpr` files

Using parmed_converter.py script to generate `.gro` and `.top` files from: - a topology file (eg an amber topology .parm7 file) - OpenMM coordinates/system information (last XML file from the production run)

```bash
python3 parmed_converter.py
```

Then;

```bash
gmx grompp -f file.mdp -c structure.gro -p topology.top -o system.tpr
```

## 2. Convert the trajectory (`.dcd` → `.xtc`)

Convert an openmm trajectory to GROMACS-compatible `.xtc` trajectory, with treating the box vectors as triclinic :

OpenMM and GROMACS use slightly different box matrix conventions.
OpenMM may contain extremely small floating-point values (e.g. `1e-16`) where GROMACS expects exact zeros..

```bash
python3 gromacs_trj.py
```

## 3. Create an index file

Generate the index file for residue that will be used for centering ;

```bash
gmx make_ndx -f structure.gro -o index.ndx
```

r RESID

#### Recommendation

For better centering, Compute the system center of mass and Select the residue closest to the center

Use:

```bash
python3 center_of_mass.py
```

## 4. Remove periodic boundary jumps and centering trajectory

```bash
gmx gromacs_nojump_centering.sh <TPR file> <TRAJ file>
```

##### Recommendation : try to not use stride ( 1/10 frames eg.) when centering ( a lots of problemes then) instead , at the end you can pass them via cpptraj to do stride , stripping of water, lipids, ions ( using strip.sh) in case of memebrane systems :

```bash
./Cpptraj_strip.sh  <TRAJ file>
```
