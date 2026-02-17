#!/usr/bin/env python3
"""
Renumber residues in PDB/CIF output files based on residue_offsets from YAML.

Boltz-2 always numbers residues starting from 1. When users provide partial
sequences with a start_index offset, this script adjusts the residue numbers
in the output structure files to match the user's original numbering.

Usage:
    renumber_structures.py <input_dir> [--yaml <yaml_file>]
    
If --yaml is not provided, looks for *.yaml in input_dir.
Only modifies files if residue_offsets are found in the YAML.
"""

import os
import re
import sys
import argparse
from glob import glob
import yaml


def load_residue_offsets(yaml_path):
    """Load residue_offsets from a YAML file."""
    with open(yaml_path, 'r') as f:
        data = yaml.safe_load(f)
    return data.get('residue_offsets', {})


def find_yaml_file(directory):
    """Find the YAML file in a directory."""
    yaml_files = glob(os.path.join(directory, '*.yaml'))
    if yaml_files:
        return yaml_files[0]
    return None


def renumber_pdb(filepath, offsets):
    """
    Renumber residues in a PDB file based on chain-specific offsets.
    
    PDB ATOM/HETATM record format (columns are 1-indexed):
    - Columns 1-6:   Record type (ATOM/HETATM)
    - Columns 7-11:  Atom serial number
    - Columns 13-16: Atom name
    - Column 17:     Alternate location indicator
    - Columns 18-20: Residue name
    - Column 22:     Chain ID
    - Columns 23-26: Residue sequence number
    - Column 27:     Code for insertion of residues
    - Columns 31-38: X coordinate
    - Columns 39-46: Y coordinate
    - Columns 47-54: Z coordinate
    """
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    modified = False
    new_lines = []
    
    for line in lines:
        if line.startswith(('ATOM  ', 'HETATM', 'TER   ')):
            # Extract chain ID (column 22, 0-indexed: 21)
            if len(line) >= 27:
                chain_id = line[21].strip()
                
                if chain_id in offsets:
                    offset = offsets[chain_id] - 1  # Convert start_index to offset
                    
                    # Extract current residue number (columns 23-26, 0-indexed: 22-25)
                    try:
                        res_num_str = line[22:26].strip()
                        res_num = int(res_num_str)
                        new_res_num = res_num + offset
                        
                        # Format: right-justified in 4 characters
                        new_res_num_str = f"{new_res_num:>4d}"
                        line = line[:22] + new_res_num_str + line[26:]
                        modified = True
                    except ValueError:
                        pass  # Skip lines with non-numeric residue numbers
        
        new_lines.append(line)
    
    if modified:
        with open(filepath, 'w') as f:
            f.writelines(new_lines)
        print(f"Renumbered PDB: {os.path.basename(filepath)}")
    
    return modified


def renumber_cif(filepath, offsets):
    """
    Renumber residues in an mmCIF file based on chain-specific offsets.
    
    In mmCIF, the relevant fields in _atom_site are:
    - auth_asym_id: Chain ID
    - auth_seq_id: Residue sequence number (user-facing)
    - label_seq_id: Residue sequence number (internal)
    """
    with open(filepath, 'r') as f:
        content = f.read()
    
    lines = content.split('\n')
    modified = False
    new_lines = []
    
    in_atom_site = False
    atom_site_columns = []
    auth_asym_idx = None
    auth_seq_idx = None
    label_seq_idx = None
    
    for line in lines:
        stripped = line.strip()
        
        # Detect start of _atom_site loop
        if stripped == 'loop_':
            in_atom_site = False
            atom_site_columns = []
            auth_asym_idx = None
            auth_seq_idx = None
            label_seq_idx = None
            new_lines.append(line)
            continue
        
        if stripped.startswith('_atom_site.'):
            in_atom_site = True
            col_name = stripped.split('.')[1].strip()
            atom_site_columns.append(col_name)
            
            if col_name == 'auth_asym_id':
                auth_asym_idx = len(atom_site_columns) - 1
            elif col_name == 'auth_seq_id':
                auth_seq_idx = len(atom_site_columns) - 1
            elif col_name == 'label_seq_id':
                label_seq_idx = len(atom_site_columns) - 1
            
            new_lines.append(line)
            continue
        
        # Process data lines within _atom_site
        if in_atom_site and auth_asym_idx is not None and stripped and not stripped.startswith('_') and not stripped.startswith('#') and not stripped == 'loop_':
            if stripped.startswith('ATOM') or stripped.startswith('HETATM'):
                parts = stripped.split()
                
                if len(parts) > max(filter(lambda x: x is not None, [auth_asym_idx, auth_seq_idx, label_seq_idx])):
                    chain_id = parts[auth_asym_idx]
                    
                    if chain_id in offsets:
                        offset = offsets[chain_id] - 1
                        
                        # Renumber auth_seq_id
                        if auth_seq_idx is not None:
                            try:
                                res_num = int(parts[auth_seq_idx])
                                parts[auth_seq_idx] = str(res_num + offset)
                                modified = True
                            except ValueError:
                                pass
                        
                        # Also renumber label_seq_id
                        if label_seq_idx is not None:
                            try:
                                res_num = int(parts[label_seq_idx])
                                parts[label_seq_idx] = str(res_num + offset)
                                modified = True
                            except ValueError:
                                pass
                    
                    new_lines.append(' '.join(parts))
                    continue
            elif not stripped:
                in_atom_site = False
        
        new_lines.append(line)
    
    if modified:
        with open(filepath, 'w') as f:
            f.write('\n'.join(new_lines))
        print(f"Renumbered CIF: {os.path.basename(filepath)}")
    
    return modified


def main():
    parser = argparse.ArgumentParser(
        description='Renumber residues in PDB/CIF files based on residue_offsets from YAML')
    parser.add_argument('input_dir', help='Directory containing structure files and YAML')
    parser.add_argument('--yaml', help='Path to YAML file with residue_offsets (auto-detected if not provided)')
    args = parser.parse_args()
    
    # Find YAML file
    yaml_path = args.yaml if args.yaml else find_yaml_file(args.input_dir)
    if not yaml_path:
        print("No YAML file found - skipping renumbering")
        sys.exit(0)
    
    # Load offsets
    offsets = load_residue_offsets(yaml_path)
    if not offsets:
        print("No residue_offsets found in YAML - skipping renumbering")
        sys.exit(0)
    
    print(f"Residue offsets: {offsets}")
    
    # Process PDB files
    pdb_files = glob(os.path.join(args.input_dir, '*.pdb'))
    for pdb_file in pdb_files:
        renumber_pdb(pdb_file, offsets)
    
    # Process CIF files
    cif_files = glob(os.path.join(args.input_dir, '*.cif'))
    for cif_file in cif_files:
        # Skip template CIF files
        if os.path.basename(cif_file).isupper():
            continue
        renumber_cif(cif_file, offsets)
    
    print("Renumbering complete")


if __name__ == '__main__':
    main()
