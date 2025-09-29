#!/usr/bin/env python3
import yaml
import re
import json
import sys
import ast
import os
from typing import List, Dict, Optional

class FlowList(list):
    pass

def flow_list_representer(dumper, data):
    return dumper.represent_sequence('tag:yaml.org,2002:seq', data, flow_style=True)
yaml.add_representer(FlowList, flow_list_representer)

PROTEIN_ALPHABET = set('ACDEFGHIKLMNPQRSTVWYUBXZ')
DNA_ALPHABET     = set('ACGTN')
RNA_ALPHABET     = set('ACGUN')

###############################
# Sequence/Ligand Type Detection
###############################
def autodetect_seq_type(seq: str) -> Optional[str]:
    s = seq.upper()
    # Remove gap, whitespace, ambiguous chars
    s = re.sub(r'[^A-Z]', '', s)
    if not s: return None
    if set(s) <= DNA_ALPHABET: return 'dna'
    if set(s) <= RNA_ALPHABET: return 'rna'
    # Heuristic: is valid protein? (don't allow only "ACGTUN" or only DNA/RNA)
    if set(s) <= PROTEIN_ALPHABET: return 'protein'
    # Most likely protein else
    return 'protein'

def validate_explicit_sequence(i, params):
    seq_key = f'seq{i}'
    seq = params.get(f'{seq_key}_seq')
    chain_id = params.get(f'{seq_key}_chain_id')
    seq_type = params.get(f'{seq_key}_type')
    # Requiredness
    if seq and not chain_id:
        sys.stderr.write(f"Error: {seq_key}_seq is set but chain_id is missing.")
        sys.exit(1)
    if chain_id and not seq:
        sys.stderr.write(f"Error: {seq_key}_chain_id is set but seq is missing.")
        sys.exit(1)
    if not seq and not chain_id:
        return None
    if not seq_type or seq_type not in ("protein", "dna", "rna"):
        sys.stderr.write(f"Error: {seq_key}_type must be specified as 'protein', 'dna', or 'rna'.")
        sys.exit(1)

    # Validate character set
    s = seq.upper()
    if seq_type == 'protein':
        valid = PROTEIN_ALPHABET
    elif seq_type == 'dna':
        valid = DNA_ALPHABET
    elif seq_type == 'rna':
        valid = RNA_ALPHABET
    else:
        sys.stderr.write(f"Error: Unknown sequence type '{seq_type}' for {seq_key}.")
        sys.exit(1)
    invalid = sorted(set(s) - valid)
    if invalid:
        sys.stderr.write(f"Error: Invalid character(s) {invalid} in {seq_type} sequence for {seq_key}.")
        sys.exit(1)
    return seq_type

def validate_no_duplicate_chain_ids_within_field(params: dict):
    for kind in ['seq', 'lig']:
        for i in range(1, 5):
            key = f'{kind}{i}_chain_id'
            cids = params.get(key)
            if cids:
                cid_list = [x.strip() for x in cids.split(',') if x.strip()]
                if len(cid_list) != len(set(cid_list)):
                    sys.stderr.write(f"Error: Duplicate chain IDs within field {key}: {cid_list}")
                    sys.exit(1)
    # Extra
    for field, block_pat in [('sequences_extra', r"\[([A-Za-z0-9_,]+)\]"), ('ligands_extra', r"\[([A-Za-z0-9_,]+)\]")]:
        extra = params.get(field)
        if extra:
            seen = set()
            for m in re.finditer(block_pat, extra):
                ids = [i.strip() for i in m.group(1).split(',')]
                for id_ in ids:
                    if id_ in seen:
                        sys.stderr.write(f"Error: Duplicate chain ID '{id_}' found in {field}.")
                        sys.exit(1)
                    seen.add(id_)

def validate_unique_chain_ids(params: Dict):
    all_chain_ids = []
    for i in range(1, 5):
        ids = params.get(f'seq{i}_chain_id')
        if ids:
            all_chain_ids.extend(x.strip() for x in ids.split(',') if x.strip())
    for i in range(1, 5):
        ids = params.get(f'lig{i}_chain_id')
        if ids:
            all_chain_ids.extend(x.strip() for x in ids.split(',') if x.strip())
    # Extra sequences
    seqe = params.get('sequences_extra')
    if seqe:
        blocks = re.findall(r"\[([A-Za-z0-9_,]+)\]", seqe)
        for ids in blocks:
            all_chain_ids.extend(x.strip() for x in ids.split(','))
    lige = params.get('ligands_extra')
    if lige:
        blocks = re.findall(r"\[([A-Za-z0-9_,]+)\]", lige)
        for ids in blocks:
            all_chain_ids.extend(x.strip() for x in ids.split(','))
    seen = set()
    duplicates = set()
    for cid in all_chain_ids:
        if cid in seen:
            duplicates.add(cid)
        else:
            seen.add(cid)
    if duplicates:
        sys.stderr.write(f"Error: Duplicate chain ID(s) found across entries: {', '.join(sorted(duplicates))}")
        sys.exit(1)

def validate_affinity_ligand_id(params: dict):
    defined_ligand_ids = set()
    for i in range(1, 5):
        chain_id = params.get(f'lig{i}_chain_id')
        if chain_id:
            for cid in [x.strip() for x in chain_id.split(',') if x.strip()]:
                defined_ligand_ids.add(cid)
    # Extra ligands
    ligextr = params.get('ligands_extra')
    if ligextr:
        for ids in re.findall(r"\[([A-Za-z0-9_,]+)\]", ligextr):
            for cid in ids.split(','):
                defined_ligand_ids.add(cid.strip())
    affinity_id = params.get('affinity_ligand_id')
    if affinity_id and affinity_id not in defined_ligand_ids:
        sys.stderr.write(f"Error: affinity_ligand_id '{affinity_id}' is not among defined ligand chain IDs: {sorted(defined_ligand_ids)}")
        sys.exit(1)

#############################
# Parsing helpers for extras
#############################
def parse_sequences_extra(seqstr: Optional[str]) -> List[Dict]:
    if not seqstr: return []
    seqs = []
    # Block: [IDs][SEQUENCE]
    pat = r"\[([A-Za-z0-9_,]+)\]\[([A-Za-z*.,\-\/\s]+)\]"
    for ids, sequence in re.findall(pat, seqstr):
        ids_list = [x.strip() for x in ids.split(',') if x.strip()]
        seq_type = autodetect_seq_type(sequence)
        if not seq_type:
            sys.stderr.write(f"Error auto-detecting sequence type for extra sequence '{sequence}'.")
            sys.exit(1)
        s = sequence.upper()
        valid = PROTEIN_ALPHABET if seq_type == "protein" else DNA_ALPHABET if seq_type == "dna" else RNA_ALPHABET
        invalid = sorted(set(s) - valid)
        if invalid:
            sys.stderr.write(f"Error: Invalid character(s) {invalid} in {seq_type} sequence for extra chains {ids}.")
            sys.exit(1)
        entry = {
            seq_type: {
                "id": ids_list[0] if len(ids_list) == 1 else FlowList(ids_list),
                "sequence": sequence,
                "msa": True if seq_type == "protein" else "empty"
            }
        }
        seqs.append(entry)
    return seqs

def parse_ligands_extra(ligstr: Optional[str]) -> List[Dict]:
    if not ligstr: return []
    ligs = []
    # Block: [IDs][CCD][SMILES]
    pat = r"\[([A-Za-z0-9_,]+)\]\[([A-Za-z0-9\-+]*)\]\[([^\[\]]*)\]"
    for ids, ccd, smiles in re.findall(pat, ligstr):
        ids_list = [x.strip() for x in ids.split(',') if x.strip()]
        if not (ccd or smiles):
            sys.stderr.write(f"Error: Must supply CCD or SMILES for ligand ids {ids}.")
            sys.exit(1)
        obj = {"ligand": {"id": ids_list[0] if len(ids_list) == 1 else FlowList(ids_list)}}
        # SMILES takes precedence
        if smiles:
            obj["ligand"]["smiles"] = smiles
        elif ccd:
            obj["ligand"]["ccd"] = ccd
        ligs.append(obj)
    return ligs


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

def parse_templates_from_params(params: Dict) -> Optional[list]:
    """
    Parses up to 2 explicit template fields (template1_path/template1_pdb_code, etc.)
    and the template_extra block string (which now contains only PDB codes).
    Returns a list of YAML template dicts, or None if nothing to set.
    """
    templates = []

    # Up to 2 explicit templates; prefer pdb_code, then path
    for i in range(1, 3):
        pdb_code = params.get(f'template{i}_pdb_code')
        path = params.get(f'template{i}_path')
        if pdb_code:
            cif_filename = f"{pdb_code.upper()}.cif"
        elif path:
            cif_filename = os.path.basename(path)
        else:
            continue

        template_entry = {'cif': cif_filename}
        # Handle the rest: prediction_id, template_id, threshold/force
        pred_ids = params.get(f'template{i}_prediction_id')
        if pred_ids:
            ids = [x.strip() for x in pred_ids.split(',') if x.strip()]
            template_entry['chain_id'] = ids[0] if len(ids) == 1 else FlowList(ids)
        templ_ids = params.get(f'template{i}_template_id')
        if templ_ids:
            ids = [x.strip() for x in templ_ids.split(',') if x.strip()]
            template_entry['template_id'] = ids[0] if len(ids) == 1 else FlowList(ids)
        thresh = params.get(f'template{i}_threshold')
        if thresh is not None:
            try:
                template_entry['force'] = True
                template_entry['threshold'] = float(thresh)
            except Exception:
                sys.stderr.write(f"Warning: could not parse threshold for template{i}: {thresh}")
        templates.append(template_entry)

    # template_extra: expect only PDB code (not a file path) in first field
    extra = params.get('template_extra')
    if extra:
        # Block format: [PDB][PREDICTION_IDS][TEMPLATE_IDS][THRESHOLD]
        block_pattern = r'(\[[^\]]+\](?:\[[^\]]*\]){0,3})'  # up to 4 fields
        blocks = re.findall(block_pattern, extra)
        for block in blocks:
            fields = re.findall(r'\[([^\]]*)\]', block)
            if not fields or not fields[0].strip():
                continue  # must have a PDB code
            pdb_code = fields[0].strip()
            # Only allow 4-character pdb codes (defensive)
            if not re.match(r'^[0-9A-Za-z]{4}$', pdb_code):
                sys.stderr.write(f"Warning: Unexpected template_extra code: '{pdb_code}'")
                continue
            entry = {'cif': f"{pdb_code.upper()}.cif"}
            if len(fields) > 1 and fields[1].strip():
                ids = [x.strip() for x in fields[1].split(',') if x.strip()]
                entry['chain_id'] = ids[0] if len(ids) == 1 else FlowList(ids)
            if len(fields) > 2 and fields[2].strip():
                ids = [x.strip() for x in fields[2].split(',') if x.strip()]
                entry['template_id'] = ids[0] if len(ids) == 1 else FlowList(ids)
            if len(fields) > 3 and fields[3].strip():
                try:
                    entry['force'] = True
                    entry['threshold'] = float(fields[3])
                except Exception:
                    sys.stderr.write(f"Warning: could not parse threshold in template_extra: '{fields[3]}'")
            templates.append(entry)

    return templates if templates else None

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

    # Explicit sequences
    for i in range(1, 5):
        seq_key = f'seq{i}'
        seq = params.get(f'{seq_key}_seq')
        chain_id = params.get(f'{seq_key}_chain_id')
        seq_type = params.get(f'{seq_key}_type')  # <- explicit
        if not seq or not chain_id or not seq_type:
            continue
        ids = [x.strip() for x in chain_id.split(',')]
        id_value = ids[0] if len(ids) == 1 else FlowList(ids)
        seq_entry = {
            seq_type: {
                'id': id_value,
                'sequence': seq
            }
        }
        mods = parse_modifications(params.get(f'{seq_key}_mod'))
        if mods:
            seq_entry[seq_type]['modifications'] = mods
        # MSA only for protein
        if seq_type == "protein":
            if params.get(f'{seq_key}_genmsa') in [None, True]:
                pass
            else:
                seq_entry[seq_type]['msa'] = "empty"
        yaml_data["sequences"].append(seq_entry)

    # Extra sequences
    seqs_extra = parse_sequences_extra(params.get('sequences_extra'))
    yaml_data["sequences"].extend(seqs_extra)

    # Explicit ligands
    for i in range(1, 5):
        lig_key = f'lig{i}'
        chain_id = params.get(f'{lig_key}_chain_id')
        ccd = params.get(f'{lig_key}_ccd')
        smiles = params.get(f'{lig_key}_smiles')
        if not chain_id:
            continue
        ids = [x.strip() for x in chain_id.split(',')]
        ligand = {'ligand': {"id": ids[0] if len(ids) == 1 else FlowList(ids)}}
        if smiles:
            ligand['ligand']['smiles'] = smiles
        elif ccd:
            ligand['ligand']['ccd'] = ccd
        else:
            sys.stderr.write(f"Error: Ligand chain_id '{chain_id}' provided but neither CCD nor SMILES given for {lig_key}.")
            sys.exit(1)
        yaml_data['sequences'].append(ligand)

    # Extra ligands
    ligs_extra = parse_ligands_extra(params.get('ligands_extra'))
    yaml_data['sequences'].extend(ligs_extra)

    # Constraints, templates, properties
    if constraints := build_constraints(params):
        yaml_data['constraints'] = constraints
    templates = parse_templates_from_params(params)
    if templates:
        yaml_data['templates'] = templates
    if properties := build_properties(params):
        yaml_data['properties'] = properties

    return yaml_data

#############################
if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.stderr.write("Error: Missing parameters")
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

        # If affinity_ligand_id is not provided, default to lig1_chain_id if available
        if not cleaned_params.get("affinity_ligand_id"):
            lig1 = cleaned_params.get("lig1_chain_id")
            if lig1:
                # Check if ligand has multiple chain IDs (comma-separated)
                chain_ids = [x.strip() for x in lig1.split(",") if x.strip()]
                if len(chain_ids) > 1:
                    sys.stderr.write("Error: Cannot compute affinity for a ligand that has multiple copies!")
                    sys.exit(1)
                # Only use the chain id from lig1_chain_id
                cleaned_params["affinity_ligand_id"] = chain_ids[0]

        # Validation
        validate_no_duplicate_chain_ids_within_field(cleaned_params)
        validate_unique_chain_ids(cleaned_params)
        validate_affinity_ligand_id(cleaned_params)
        # Per-sequence validation
        for i in range(1, 5):
            validate_explicit_sequence(i, cleaned_params)

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
        sys.stderr.write(f"Error decoding JSON: {str(e)}")
        sys.exit(1)