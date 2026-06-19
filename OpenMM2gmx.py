"""
OpenMM → GROMACS Trajectory Conversion and Centering Workflow
=============================================================

Author: Alaa REGUEI

Overview
--------
This script provides a complete workflow to convert and process molecular dynamics
trajectories generated with OpenMM into a format compatible with GROMACS analysis tools.

It is designed to handle systems with periodic boundary conditions, including triclinic
boxes, and has been tested on complex biomolecular systems such as membrane proteins.

The workflow includes:
    1. Generation of GROMACS-compatible topology and coordinate files (.gro, .top, .tpr)
       from OpenMM inputs (.parm7 and system XML).
    2. Conversion of OpenMM .dcd trajectories into GROMACS-readable .xtc format.
       Special care is taken to correctly handle differences in box definitions between
       OpenMM and GROMACS.
    3. Post-processing of trajectories to remove periodic boundary jumps.
    4. Centering of trajectories around a user-defined selection (e.g., protein COM or
       specific residues).

Features
--------
- Multi-trajectory support (concatenation of multiple .dcd files)
- Flexible centering options:
    - Center of mass (COM)
    - Specific residues or residue ranges
- Trajectory stride control
- Multiple output modes:
    - Full system trajectory (protein + solvent/non-polymers)
    - Protein-only trajectory
    - Both outputs
- Handling of periodic boundary conditions

Dependencies
------------
- GROMACS (must be installed and accessible via `gmx` command)
- Python dependencies listed in requirements.txt

"""

import parmed as pmd
import subprocess
import MDAnalysis as mda
from MDAnalysis.coordinates.XTC import XTCWriter
import warnings
warnings.filterwarnings('ignore')
import numpy as np
import os
import sys
import argparse
import logging
from tqdm import tqdm
from openmm.app import PDBFile, PDBxFile

# Logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
logger.addHandler(logging.StreamHandler(sys.stdout))


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--top", dest='topfile', required=True, help='Input topology file : eg. parm7 or pdb ',)
    parser.add_argument("--xml", dest='xmlfile', required=True, help='Input OpenMM coordinates/system information XML file from the production run ',)
    parser.add_argument("--traj", dest='trajfile', nargs="+", required=True, help='Input OpenMM trajectory files : eg. traj.dcd or traj_part1.dcd traj_part2.dcd ...')
    parser.add_argument("--mdp", dest='mdpfile', required=True, help='Input GROMACS mdp file',)
    parser.add_argument("--center_res", dest='center_res', nargs="+", required=True, help='Residue(s) used to center the simulation. Use "COM" to center on the protein center of mass, or specify one or more residues as ranges in this form "resid1-resid3" or "resid1" if single residues. Example: --center_res COM or --center_res 10-50 or --center_res 10 20 30',)
    parser.add_argument("--sim_name", dest='sim_name', required=True, help='Name of the simulation_replica if needed',)
    parser.add_argument("--stride", dest='stride', required=False, help='Stride for trajectory analysis , not recommended if the trajectory is partitioned', type=int , default=1)
    parser.add_argument("--save_mode", dest="save_mode", type=int, choices=[0, 1, 2], default=0, help=("0 = full corrected system trajectory (Protein + non-polymer stuff)," \
        "1 = protein-only corrected trajectory, "\
        "2 = both trajectories"),)   
    parser.add_argument("--output_dir", dest='output_dir', default="MD_output", help='Output directory name',)
    return parser.parse_args()

def generate_gro_top_files(top_file, xml_file, sim_name, output_dir):
    logger.info("- Generating .gro and .top files")
    if top_file.endswith(".parm7"):
        system = pmd.load_file(top_file, xyz=xml_file)
    elif top_file.endswith(".pdb"):
        system = pmd.openmm.load_topology(PDBFile(top_file).topology, xyz=xml_file)
    
    # Output files
    GRO_OUTPUT = f"{output_dir}/{sim_name}.gro"
    TOP_OUTPUT = f"{output_dir}/{sim_name}.top"
    
    logger.info ("- Detected box vectors: %s", system.box)
    system.save(GRO_OUTPUT, format="gro", overwrite=True)
    system.save(TOP_OUTPUT, format="gromacs", overwrite=True)
    with open(f"{output_dir}/log_files/{sim_name}_generation_gro_top.log", "w") as log_file:
        log_file.write(f"Detected box vectors:\n{system.box}\n")
        log_file.write(f"Generated files: {GRO_OUTPUT}, {TOP_OUTPUT}\n")
        
    logger.info(f"- Conversion completed successfully, Generated files: {GRO_OUTPUT}, {TOP_OUTPUT}")
    
    return GRO_OUTPUT, TOP_OUTPUT

def generate_tpr(GRO_OUTPUT, TOP_OUTPUT, mdp_file, sim_name, output_dir):
    
    logger.info("- Generating TPR file using GROMACS")
    TPR_OUTPUT = f"{output_dir}/{sim_name}.tpr"
    cmd = [
        "gmx", "grompp",
        "-f", mdp_file,
        "-c", GRO_OUTPUT,
        "-p", TOP_OUTPUT,
        "-o", TPR_OUTPUT,
        "-maxwarn", "1"
    ]
    result = subprocess.run(cmd, text=True, capture_output=True)
    with open(f"{output_dir}/log_files/{sim_name}_generation_tpr.log", "w") as log_file:
        log_file.write(result.stdout)
        log_file.write(result.stderr)
        log_file.write(f"TPR file generated: {TPR_OUTPUT}")
    logger.info(f"- TPR file generated: {TPR_OUTPUT}")
    return TPR_OUTPUT
    

def convert_traj(top_file, traj_files, sim_name, stride=1, output_dir="MD_output"):

    logger.info("- Converting trajectory files to XTC format")
    if top_file.endswith(".parm7"):
        u = mda.Universe(top_file, traj_files)
    elif top_file.endswith(".pdb"):
        u = mda.Universe(PDBFile(top_file).topology, traj_files)
    with XTCWriter(f"{output_dir}/{sim_name}_trj.xtc", n_atoms=u.atoms.n_atoms) as writer:
        for ts in tqdm(u.trajectory[::stride]):
            box = ts.dimensions.copy()
            box[np.abs(box) < 1e-2] = 0.0
            ts.dimensions = box
            writer.write(u.atoms)
    OUTPUT_XTC = f"{output_dir}/{sim_name}_trj.xtc"
    logger.info(f"- Trajectory conversion completed: {sim_name}_trj.xtc")
    return OUTPUT_XTC

def compute_center_of_mass(GRO_OUTPUT):
    
    u = mda.Universe(GRO_OUTPUT)
    protein = u.select_atoms("protein")
    com_protein = protein.center_of_mass()
    logger.info(f"- Center of mass: {com_protein}")

    a = None
    min_dist = float("inf")
    for residu in protein.residues:
        com_residu = residu.atoms.center_of_mass()
        distance = np.linalg.norm(com_residu - com_protein)

        if distance < min_dist:
            min_dist = distance
            closest_residue = residu
    logger.info(
        f"- The residue closest to the center of mass is : {closest_residue.resname} (Resid: {closest_residue.resid})")
    return int(closest_residue.resid)


def make_index_file(GRO_OUTPUT, center_res, sim_name, output_dir):
    logger.info("- Generating index file ")
    cmd = ["gmx", "make_ndx", "-f", GRO_OUTPUT, "-o", f"{sim_name}_index.ndx"]

    if isinstance(center_res, list):
        center_res = center_res[0]

    if center_res == "COM":

        COM_res = compute_center_of_mass(GRO_OUTPUT)
        ndx_commands = f"""
        r {COM_res}
        name  COM
        q
        """

    elif center_res.find("-") != -1:
        ndx_commands = f"""
        r {center_res}
        name {center_res}
        q
        """

    else:
        logger.error("Invalid center_res format, please refert to documentation for correct usage")
        return

    result = subprocess.run(
        cmd,
        input=ndx_commands,
        text=True,
        capture_output=True
    )


    lines = result.stdout.splitlines()
    group_ids = []
    for line in lines:
        parts = line.strip().split()
        if len(parts) > 0 and parts[0].isdigit():
            group_ids.append(int(parts[0]))

    new_group_id = max(group_ids)+1

    # Extract the group index for the "Protein" group, used for centering in the next step if needed
    for line in lines:
        if "Protein             " in line.strip():
            protein_index = int(line.strip()[0])
            break

    with open(f"{output_dir}/log_files/{sim_name}_generation_index.log", "w") as log_file:
        log_file.write(result.stdout)
        log_file.write(result.stderr)
        log_file.write("Index file generated: index.ndx")

    logger.info(f"- Index file generated: {sim_name}_index.ndx with group {center_res} (group id: {new_group_id})")
    return new_group_id, f"{sim_name}_index.ndx", protein_index



def traj_correction(GRO_OUTPUT, TPR_OUTPUT, OUTPUT_XTC, INDEX_OUTPUT, new_group_id, sim_name, output_dir, save_mode, protein_index):
    logger.info("- Correcting trajectory for PBC and centering")

    nojump_xtc = f"{output_dir}/{sim_name}_nojump.xtc"
    centered_xtc_full_complex = f"{output_dir}/{sim_name}_trj_no_jump_centered_full_complex.xtc"
    centered_xtc_only_protein = f"{output_dir}/{sim_name}_trj_no_jump_centered_only_protein.xtc"

    log_path = f"{output_dir}/log_files/{sim_name}_traj_no_jump_centered.log"

    with open(log_path, "w") as log_file:
        # Step 1: remove jumps
        cmd1 = ["gmx", "trjconv","-s", TPR_OUTPUT,"-f", OUTPUT_XTC,"-o", nojump_xtc,"-pbc", "nojump"]

        result1 = subprocess.run( cmd1, input="0\n", text=True, capture_output=True)
        log_file.write("=== Step 1: -pbc nojump ===\n")
        log_file.write(result1.stdout)
        log_file.write(result1.stderr)
        log_file.write("\n")

        if result1.returncode != 0:
            logger.error("Step 1 failed during trajectory nojump correction.")
            raise RuntimeError("gmx trjconv failed at step 1 (-pbc nojump).")

        # Step 2: center and make compact
        
        if save_mode == 0:
            centered_xtc = centered_xtc_full_complex
            cmd2 = ["gmx", "trjconv","-s", TPR_OUTPUT,"-f", nojump_xtc,"-n", INDEX_OUTPUT,"-o", centered_xtc,"-pbc", "mol","-center","-ur", "compact"]
            result2 = subprocess.run( cmd2, input=f"{new_group_id}\n0\n", text=True, capture_output=True)        
        elif save_mode == 1:
            centered_xtc = centered_xtc_only_protein
            cmd2 = ["gmx", "trjconv","-s", TPR_OUTPUT,"-f", nojump_xtc,"-n", INDEX_OUTPUT,"-o", centered_xtc,"-pbc", "mol","-center","-ur", "compact"]
            result2 = subprocess.run( cmd2, input=f"{new_group_id}\n{protein_index}\n", text=True, capture_output=True)        
            
            cmd3 = ["gmx", "trjconv","-s", TPR_OUTPUT,"-f", GRO_OUTPUT,"-n", INDEX_OUTPUT,"-o", f"{output_dir}/{sim_name}.gro"]
            result3 = subprocess.run( cmd3, input=f"{protein_index}\n", text=True, capture_output=True)        

            
        elif save_mode == 2:
            centered_xtc = centered_xtc_full_complex
            cmd2 = ["gmx", "trjconv","-s", TPR_OUTPUT,"-f", nojump_xtc,"-n", INDEX_OUTPUT,"-o", centered_xtc,"-pbc", "mol","-center","-ur", "compact"]
            result2 = subprocess.run( cmd2, input=f"{new_group_id}\n0\n", text=True, capture_output=True)        
            
            centered_xtc = centered_xtc_only_protein
            cmd2 = ["gmx", "trjconv","-s", TPR_OUTPUT,"-f", nojump_xtc,"-n", INDEX_OUTPUT,"-o", centered_xtc,"-pbc", "mol","-center","-ur", "compact"]
            result2 = subprocess.run( cmd2, input=f"{new_group_id}\n{protein_index}\n", text=True, capture_output=True) 
            cmd3 = ["gmx", "trjconv","-s", TPR_OUTPUT,"-f", GRO_OUTPUT,"-n", INDEX_OUTPUT,"-o", f"{output_dir}/{sim_name}.gro"]
            result3 = subprocess.run( cmd3, input=f"{protein_index}\n", text=True, capture_output=True)        
       

        log_file.write("=== Step 2: centering ===\n")
        log_file.write(result2.stdout)
        log_file.write(result2.stderr)
        log_file.write(f"\nCorrected trajectory generated: {centered_xtc}\n")

        if result2.returncode != 0:
            logger.error("Step 2 failed during centering.")
            raise RuntimeError("gmx trjconv failed at step 2 (-center).")

    logger.info(f"- Trajectory correction completed ! ")
    return centered_xtc


if __name__ == "__main__":

    args = parse_args()
    logger.info("- Read input files")
    top_file = args.topfile
    xml_file = args.xmlfile
    traj_files = args.trajfile
    mdp_file = args.mdpfile
    sim_name = args.sim_name
    center_res = args.center_res
    stride = args.stride
    save_mode = args.save_mode
    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(f"{output_dir}/log_files", exist_ok=True)
    
    
    GRO_OUTPUT, TOP_OUTPUT = generate_gro_top_files(top_file, xml_file, sim_name, output_dir)
    TPR_OUTPUT = generate_tpr(GRO_OUTPUT, TOP_OUTPUT, mdp_file, sim_name, output_dir)

    OUTPUT_XTC = convert_traj(top_file, traj_files, sim_name, stride, output_dir)
    new_group_id, INDEX_OUTPUT, protein_index = make_index_file(GRO_OUTPUT, center_res, sim_name, output_dir)
    traj_correction(GRO_OUTPUT, TPR_OUTPUT, OUTPUT_XTC, INDEX_OUTPUT, new_group_id, sim_name, output_dir, save_mode, protein_index)
