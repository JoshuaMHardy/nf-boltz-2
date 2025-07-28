#!/usr/bin/env python
import yaml
import re
import json
import sys
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

def build_constraints(params: Dict) -> Optional[List[Dict]]:
    constraints = []
    
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
    
    if params.get('pocket_ligand_chain_id') and params.get('pocket_contacts'):
        contact_pairs = []
        for contact in parse_constraint_atoms(params['pocket_contacts']):
            if len(contact) >= 2:
                contact_pairs.append(FlowList(contact[:2]))  # Each contact as FlowList
        if contact_pairs:
            # Wrap the entire contacts list in FlowList for single-line formatting
            constraints.append({
                'pocket': {
                    'binder': params['pocket_ligand_chain_id'].strip(),
                    'contacts': FlowList(contact_pairs)
                }
            })
    
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
        if not params.get(f'{lig_key}_type'):
            continue
            
        ids = [x.strip() for x in params[f'{lig_key}_chain_id'].split(',')] if params[f'{lig_key}_chain_id'] else []
        
        # Use FlowList for multiple IDs, single string for single ID
        id_value = ids[0] if len(ids) == 1 else FlowList(ids)
        
        ligand = {
            'ligand': {
                'id': id_value
            }
        }
                
        lig_type = params[f'{lig_key}_type'].lower()
        if lig_type == 'ccd':
            ligand['ligand']['ccd'] = params[f'{lig_key}_chem']
        elif lig_type == 'smiles':
            ligand['ligand']['smiles'] = params[f'{lig_key}_chem']
        
        yaml_data['sequences'].append(ligand)
    
    if constraints := build_constraints(params):
        yaml_data['constraints'] = constraints
    
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
