#!/usr/bin/env python
import yaml
import re
import json
import sys
import ast
from typing import List, Dict, Optional

class FlowList(list):
    """Custom list class to mark lists that should use flow style"""
    pass

def flow_list_representer(dumper, data):
    """Custom representer for FlowList to force flow style"""
    return dumper.represent_sequence('tag:yaml.org,2002:seq', data, flow_style=True)

# Register the custom representer
yaml.add_representer(FlowList, flow_list_representer)

def parse_modifications(mod_str: Optional[str]) -> List[Dict]:
    if not mod_str:
        return []
    modifications = []
    pattern = r'\[(\d+),([A-Z]+)\]'
    matches = re.findall(pattern, mod_str)
    for position, ccd in matches:
        modifications.append({'position': int(position), 'ccd': ccd})
    return modifications

def parse_constraint_atoms(constraint_str: str) -> List[List[str]]:
    atoms = []
    matches = re.findall(r'\[([^]]+)\]', constraint_str)
    for match in matches:
        parts = [p.strip() for p in match.split(',')]
        if len(parts) >= 2:
            parts[1] = int(parts[1])
        atoms.append(FlowList(parts))  # Use FlowList for atom coordinates
    return atoms

def parse_pocket_constraint(pocket_str: str):
    """
    Parses a pocket constraint string into components.
    Format: '[B][[A,53],[A,57],[A,61]] [C][[A,101],[A,130]][5][true]'
    Returns a list of dicts, one per pocket constraint block.
    """
    constraints = []
    if not pocket_str:
        return constraints
    # Split multiple pockets on spaces between top-level '[...]' blocks
    pattern = r'(\[[^\]]+\]\[\[.*?\]\](?:\[\d+\])?(?:\[\w+\])?)'
    matches = re.findall(pattern, pocket_str)
    for match in matches:
        # Extract binder
        binder = re.findall(r'^\[([^\]]+)\]', match)
        # Extract contacts
        contacts = re.findall(r'\[\[([^\]]+)\]\]', match)
        # Extract max_distance (optional)
        dist_match = re.search(r'\[(\d+(?:\.\d+)?)\](?![^\[]*\])', match)
        max_distance = float(dist_match.group(1)) if dist_match else None
        # Extract force (optional)
        force_match = re.search(r'\[(true|false)\](?![^\[]*\])', match, re.IGNORECASE)
        force = force_match.group(1).lower() == "true" if force_match else None

        if binder and contacts:
            contacts_list = [ast.literal_eval(f'[{c.strip()}]') for c in contacts[0].split('],[')]
            constraint = {
                'pocket': {
                    'binder': binder[0],
                    'contacts': [FlowList(x) for x in contacts_list]
                }
            }
            if max_distance is not None:
                constraint['pocket']['max_distance'] = max_distance
            if force is not None:
                constraint['pocket']['force'] = force
            constraints.append(constraint)
    return constraints

def parse_contact_constraint(contact_str: str):
    """
    Parses a contact constraint string to extract atom pairs, optional distance, and force.
    Format: '[A,53][A,57][5][true] [B,101][A,130][4.5][false]'
    Returns a list of dicts, one per contact block.
    """
    constraints = []
    if not contact_str:
        return constraints
    # Matches blocks like: [A,53][A,57][5][true]
    pattern = r'(\[[^\]]+\]\[[^\]]+\](?:\[\d+(?:\.\d+)?\])?(?:\[(?:true|false)\])?)'
    matches = re.findall(pattern, contact_str, re.IGNORECASE)
    for match in matches:
        atoms = re.findall(r'\[([^\]]+)\]', match)
        atom1 = [x.strip() for x in atoms[0].split(',')]
        atom2 = [x.strip() for x in atoms[1].split(',')]
        for i in [1, 1]:
            if len(atom1) > 1:
                atom1[1] = int(atom1[1])
            if len(atom2) > 1:
                atom2[1] = int(atom2[1])

        # Optional max_distance
        dist_match = re.search(r'\[(\d+(?:\.\d+)?)\](?![^\[]*\])', match)
        max_distance = float(dist_match.group(1)) if dist_match else None

        # Optional force flag
        force_match = re.search(r'\[(true|false)\](?![^\[]*\])', match, re.IGNORECASE)
        force = force_match.group(1).lower() == "true" if force_match else None

        constraint = {
            'contact': {
                'atom1': FlowList(atom1),
                'atom2': FlowList(atom2)
            }
        }
        if max_distance is not None:
            constraint['contact']['max_distance'] = max_distance
        if force is not None:
            constraint['contact']['force'] = force
        constraints.append(constraint)
    return constraints

def parse_templates(templates_str: Optional[str]) -> Optional[List[Dict]]:
    """
    Parse templates string as specified for Boltz-2
    """
    if not templates_str:
        return None
    # Find all '[...]' fields
    template_blocks = []
    # Split template blocks: a block starts with [PATH], there can be multiple in the string
    block_pattern = r'(\[[^\]]+\](?:\[[^\]]*\]){0,4})'  # up to 5 fields
    blocks = re.findall(block_pattern, templates_str)
    for block in blocks:
        # Extract all [xxx] fields
        fields = re.findall(r'\[([^\]]*)\]', block)
        if not fields:
            continue
        entry = {
            'cif': fields[0]  # Always exists
        }
        # chain_id(s)
        if len(fields) > 1 and fields[1].strip():
            ids = [x.strip() for x in fields[1].split(',') if x.strip()]
            if len(ids) == 1:
                entry['chain_id'] = ids[0]
            elif len(ids) > 1:
                entry['chain_id'] = FlowList(ids)
        # template_id(s)
        if len(fields) > 2 and fields[2].strip():
            ids = [x.strip() for x in fields[2].split(',') if x.strip()]
            if len(ids) == 1:
                entry['template_id'] = ids[0]
            elif len(ids) > 1:
                entry['template_id'] = FlowList(ids)
        # force (optional, default False)
        if len(fields) > 3 and fields[3].strip():
            entry['force'] = fields[3].lower() == 'true'
        # threshold (optional, float or int)
        if len(fields) > 4 and fields[4].strip():
            thresh = fields[4]
            entry['threshold'] = float(thresh) if '.' in thresh else int(thresh)
        template_blocks.append(entry)
    return template_blocks if template_blocks else None

def build_properties(params: Dict) -> Optional[list]:
    """
    Builds the 'properties' block for the YAML file.
    For now, only supports affinity binder.
    Returns a list, or None if nothing to set.
    """
    props = []
    binder_id = params.get('affinity_ligand_id')
    if binder_id:
        props.append({'affinity': {'binder': binder_id}})
    return props if props else None

def build_constraints(params: Dict) -> Optional[List[Dict]]:
    constraints = []

    # Bond constraints
    if params.get('constraints_bond'):
        for bond_pair in params['constraints_bond'].split():
            atoms = parse_constraint_atoms(bond_pair)
            if len(atoms) == 2:
                constraints.append({
                    'bond': {
                        'atom1': atoms[0],
                        'atom2': atoms[1]
                    }
                })
    
    # Pocket constraints
    if params.get('constraints_pocket'):
        constraints.extend(parse_pocket_constraint(params['constraints_pocket']))  # New parser

    # Contact constraints
    if params.get('constraints_contacts'):
        constraints.extend(parse_contact_constraint(params['constraints_contacts']))  # New parser

    return constraints if constraints else None


def build_yaml_structure(params: Dict) -> Dict:
    yaml_data = {'sequences': []}
    
    for i in range(1, 13):
        seq_key = f'seq{i}'
        if not params.get(f'{seq_key}_type'):
            continue
            
        seq_type = params[f'{seq_key}_type']
        ids = [x.strip() for x in params[f'{seq_key}_chain_id'].split(',')] if params[f'{seq_key}_chain_id'] else []
        
        # Use FlowList for multiple IDs, single string for single ID
        id_value = ids[0] if len(ids) == 1 else FlowList(ids)
        
        seq_entry = {
            seq_type: {
                'id': id_value,
                'sequence': params[f'{seq_key}_seq']
            }
        }
        
        if mods := parse_modifications(params[f'{seq_key}_mod']):
            seq_entry[seq_type]['modifications'] = mods
        
        if seq_type in ('protein', 'dna', 'rna') and not params[f'{seq_key}_genmsa']:
            seq_entry[seq_type]['msa'] = 'empty'
            
        yaml_data['sequences'].append(seq_entry)
    
    for i in range(1, 13):
        lig_key = f'lig{i}'
        chain_id = params.get(f'{lig_key}_chain_id')
        ccd = params.get(f'{lig_key}_ccd')
        smiles = params.get(f'{lig_key}_smiles')

        if not chain_id:
            continue  # No ligand defined in this slot

        ids = [x.strip() for x in chain_id.split(',')] if chain_id else []
        id_value = ids[0] if len(ids) == 1 else FlowList(ids)
        ligand = {
            'ligand': {
                'id': id_value
            }
        }

        # Precedence: SMILES > CCD
        if smiles:
            ligand['ligand']['smiles'] = smiles
        elif ccd:
            ligand['ligand']['ccd'] = ccd
        else:
            print(f"Error: Ligand chain_id '{chain_id}' provided but neither 'CCD' nor 'SMILES' given for {lig_key}.")
            sys.exit(1)

        yaml_data['sequences'].append(ligand)
    
    if constraints := build_constraints(params):
        yaml_data['constraints'] = constraints
    
    if templates := parse_templates(params.get('templates')):
        yaml_data['templates'] = templates
    
    if properties := build_properties(params):
        yaml_data['properties'] = properties

    return yaml_data

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Error: Missing parameters")
        sys.exit(1)
        
    try:
        input_params = json.loads(sys.argv[1])
        cleaned_params = {}
        for k, v in input_params.items():
            if isinstance(v, str):
                if v.lower() == "null": v = None
                elif v.lower() == "true": v = True
                elif v.lower() == "false": v = False
            cleaned_params[k] = v
        
        yaml_structure = build_yaml_structure(cleaned_params)
        print(yaml.dump(
            yaml_structure, 
            sort_keys=False, 
            width=120, 
            default_flow_style=False,
            indent=4,
            allow_unicode=True
        ))
        
    except json.JSONDecodeError as e:
        print(f"Error decoding JSON: {str(e)}")
        sys.exit(1)
