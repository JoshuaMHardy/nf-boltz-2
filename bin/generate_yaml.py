#!/usr/bin/env python3
import yaml
import re
import json
import sys
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
    # Extra sequences: first bracket field is chain IDs
    seqe = params.get('sequences_extra')
    if seqe:
        seen = set()
        for m in re.finditer(r"\[([A-Za-z0-9_,]+)\]\[", seqe):
            ids = [x.strip() for x in m.group(1).split(',') if x.strip()]
            for id_ in ids:
                if id_ in seen:
                    sys.stderr.write(f"Error: Duplicate chain ID '{id_}' found in sequences_extra.")
                    sys.exit(1)
                seen.add(id_)
    # Extra ligands: use the proper bracket-depth parser to get only chain IDs
    ligee = params.get('ligands_extra')
    if ligee:
        seen = set()
        pos = 0
        while pos < len(ligee):
            if ligee[pos] != '[':
                pos += 1
                continue
            ids_str, pos = _read_bracket_field(ligee, pos)
            if ids_str is None:
                break
            # Skip CCD and SMILES fields
            _, pos = _read_bracket_field(ligee, pos)
            _, pos = _read_bracket_field(ligee, pos)
            for id_ in [x.strip() for x in ids_str.split(',') if x.strip()]:
                if id_ in seen:
                    sys.stderr.write(f"Error: Duplicate chain ID '{id_}' found in ligands_extra.")
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
    # Extra sequences: first bracket field per block is chain IDs
    seqe = params.get('sequences_extra')
    if seqe:
        for m in re.finditer(r"\[([A-Za-z0-9_,]+)\]\[", seqe):
            all_chain_ids.extend(x.strip() for x in m.group(1).split(',') if x.strip())
    # Extra ligands: use bracket-depth parser to get only the chain IDs field
    lige = params.get('ligands_extra')
    if lige:
        pos = 0
        while pos < len(lige):
            if lige[pos] != '[':
                pos += 1
                continue
            ids_str, pos = _read_bracket_field(lige, pos)
            if ids_str is None:
                break
            _, pos = _read_bracket_field(lige, pos)  # skip CCD
            _, pos = _read_bracket_field(lige, pos)  # skip SMILES
            all_chain_ids.extend(x.strip() for x in ids_str.split(',') if x.strip())
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
    # Block: [IDs][SEQUENCE] or [IDs][SEQUENCE][START_INDEX]
    pat = r"\[([A-Za-z0-9_,]+)\]\[([A-Za-z*.,\-\/\s]+)\](?:\[(\d+)\])?"
    for match in re.finditer(pat, seqstr):
        ids = match.group(1)
        sequence = match.group(2)
        start_index_str = match.group(3)
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
        # Store start_index if provided (will be collected later)
        if start_index_str:
            entry['_start_index'] = int(start_index_str)
        seqs.append(entry)
    return seqs

def _read_bracket_field(s: str, pos: int) -> tuple:
    """
    Read a single bracket-delimited field from string s starting at pos.
    Handles nested brackets (e.g. SMILES like [Na+], [Ca+2]).
    Returns (field_content, next_pos) or (None, pos) if no bracket at pos.
    """
    if pos >= len(s) or s[pos] != '[':
        return None, pos
    depth = 0
    start = pos + 1
    i = pos
    while i < len(s):
        if s[i] == '[':
            depth += 1
        elif s[i] == ']':
            depth -= 1
            if depth == 0:
                return s[start:i], i + 1
        i += 1
    return None, pos  # Unbalanced brackets

def parse_ligands_extra(ligstr: Optional[str]) -> List[Dict]:
    """Parse ligands_extra string. Uses bracket-depth parsing so SMILES
    containing brackets (e.g. [Na+], [Ca+2], CC(=O)[O-]) are handled.
    Format: '[IDs][CCD][SMILES]' per block, comma-separated. Whitespace is ignored.
    Example: '[G,H][ZN][],[I][][[Na+]]'
    """
    if not ligstr: return []
    ligs = []
    i = 0
    while i < len(ligstr):
        # Skip whitespace between blocks
        if ligstr[i] != '[':
            i += 1
            continue
        # Read three bracket-delimited fields: [IDs][CCD][SMILES]
        ids_str, i = _read_bracket_field(ligstr, i)
        if ids_str is None:
            break
        ccd, i = _read_bracket_field(ligstr, i)
        if ccd is None:
            sys.stderr.write(f"Error: Missing CCD field for ligand ids '{ids_str}' in ligands_extra.")
            sys.exit(1)
        smiles, i = _read_bracket_field(ligstr, i)
        if smiles is None:
            sys.stderr.write(f"Error: Missing SMILES field for ligand ids '{ids_str}' in ligands_extra.")
            sys.exit(1)
        ids_list = [x.strip() for x in ids_str.split(',') if x.strip()]
        if not (ccd or smiles):
            sys.stderr.write(f"Error: Must supply CCD or SMILES for ligand ids {ids_str}.")
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

def _parse_contact_pair(pair_str: str) -> list:
    """
    Parse a single contact pair string like 'A,53' or 'A,CA' into [chain_id, residue_idx_or_atom].
    Handles both integer residue indices and string atom names.
    """
    parts = [p.strip() for p in pair_str.split(',')]
    if len(parts) >= 2:
        try:
            parts[1] = int(parts[1])
        except ValueError:
            pass  # Keep as string (atom name for non-polymers)
    return parts


def parse_pocket_constraint(pocket_str: str):
    """
    Parses a pocket constraint string into components.
    Format: '[B][[A,53],[A,57],[A,61]],[C][[A,101],[A,130]][5][true]'
    Multiple pockets separated by commas. Whitespace is ignored.
    Returns a list of dicts, one per pocket constraint block.
    """
    constraints = []
    if not pocket_str:
        return constraints
    # Split multiple pockets on spaces between top-level '[...]' blocks
    pattern = r'(\[[^\]]+\]\[\[.*?\]\](?:\[\d+(?:\.\d+)?\])?(?:\[(?:true|false)\])?)'
    matches = re.findall(pattern, pocket_str, re.IGNORECASE)
    for match in matches:
        # Extract binder
        binder = re.findall(r'^\[([^\]]+)\]', match)
        # Extract contacts block - everything between [[ and ]]
        contacts_match = re.search(r'\[\[(.*?)\]\]', match)
        # Extract max_distance (optional) - number after the contacts block
        dist_match = re.search(r'\]\]\[(\d+(?:\.\d+)?)\]', match)
        max_distance = float(dist_match.group(1)) if dist_match else None
        # Extract force (optional)
        force_match = re.search(r'\[(true|false)\]\s*$', match, re.IGNORECASE)
        force = force_match.group(1).lower() == "true" if force_match else None

        if binder and contacts_match:
            # Parse contacts: split on ],[  to get individual [chain,idx] pairs
            contacts_raw = contacts_match.group(1)
            contact_pairs = contacts_raw.split('],[')
            contacts_list = [_parse_contact_pair(p.strip('[] ')) for p in contact_pairs]
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
    Parses a contact constraint string to extract token pairs, optional distance, and force.
    Format: '[A,53][A,57],[A,101][A,130][5][true]'
    Multiple contacts separated by commas. Whitespace is ignored.
    Returns a list of dicts, one per contact block.
    Note: Boltz-2 uses 'token1'/'token2' field names for contact constraints.
    """
    constraints = []
    if not contact_str:
        return constraints
    # Matches blocks like: [A,53][A,57][5][true]
    pattern = r'(\[[^\]]+\]\[[^\]]+\](?:\[\d+(?:\.\d+)?\])?(?:\[(?:true|false)\])?)'
    matches = re.findall(pattern, contact_str, re.IGNORECASE)
    for match in matches:
        atoms = re.findall(r'\[([^\]]+)\]', match)
        token1 = _parse_contact_pair(atoms[0])
        token2 = _parse_contact_pair(atoms[1])

        # Optional max_distance
        dist_match = re.search(r'\[(\d+(?:\.\d+)?)\](?![^\[]*\])', match)
        max_distance = float(dist_match.group(1)) if dist_match else None

        # Optional force flag
        force_match = re.search(r'\[(true|false)\](?![^\[]*\])', match, re.IGNORECASE)
        force = force_match.group(1).lower() == "true" if force_match else None

        constraint = {
            'contact': {
                'token1': FlowList(token1),
                'token2': FlowList(token2)
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


def collect_residue_offsets(params: Dict, sequences: List[Dict]) -> Dict[str, int]:
    """
    Collect start_index offsets for each chain ID from explicit params and extra sequences.
    Returns a dict mapping chain_id -> start_index (only for non-default values, i.e. != 1).
    """
    offsets = {}
    # Explicit sequences (seq1-seq4)
    for i in range(1, 5):
        chain_id = params.get(f'seq{i}_chain_id')
        start_index = params.get(f'seq{i}_start_index')
        if chain_id and start_index is not None:
            start_idx = int(start_index)
            if start_idx != 1:
                # Apply to all chain IDs (comma-separated for multi-copy)
                for cid in [x.strip() for x in chain_id.split(',') if x.strip()]:
                    offsets[cid] = start_idx
    
    # Extra sequences (may have _start_index from parsing)
    for entry in sequences:
        start_idx = entry.pop('_start_index', None)
        if start_idx and start_idx != 1:
            # Find the sequence type key (protein/dna/rna)
            for key in ('protein', 'dna', 'rna'):
                if key in entry:
                    seq_id = entry[key].get('id')
                    if isinstance(seq_id, list):
                        for cid in seq_id:
                            offsets[cid] = start_idx
                    elif seq_id:
                        offsets[seq_id] = start_idx
                    break
    
    return offsets


def adjust_modification_positions(modifications: List[Dict], start_index: int) -> List[Dict]:
    """
    Adjust modification positions from user numbering to Boltz-2 numbering.
    User provides positions in their real numbering (e.g., residue 35 with start_index=30).
    Boltz-2 expects 1-based indexing relative to the provided sequence.
    boltz_pos = user_pos - start_index + 1
    """
    adjusted = []
    for mod in modifications:
        new_mod = dict(mod)
        new_mod['position'] = mod['position'] - start_index + 1
        if new_mod['position'] < 1:
            sys.stderr.write(
                f"Error: Modification position {mod['position']} is before the sequence start index {start_index}.")
            sys.exit(1)
        adjusted.append(new_mod)
    return adjusted


def adjust_constraint_residue_index(chain_id: str, res_idx: int, offsets: Dict[str, int]) -> int:
    """
    Adjust a single constraint residue index from user numbering to Boltz-2 numbering.
    Returns the adjusted 1-based index for Boltz-2.
    """
    start_index = offsets.get(chain_id, 1)
    if start_index == 1:
        return res_idx
    adjusted = res_idx - start_index + 1
    if adjusted < 1:
        sys.stderr.write(
            f"Error: Constraint residue index {res_idx} for chain {chain_id} is before the sequence start index {start_index}.")
        sys.exit(1)
    return adjusted


def adjust_constraints(constraints: List[Dict], offsets: Dict[str, int]) -> List[Dict]:
    """
    Adjust all constraint residue indices from user numbering to Boltz-2 numbering.
    Handles bond, contact, and pocket constraints.
    """
    if not offsets or not constraints:
        return constraints
    
    adjusted = []
    for constraint in constraints:
        new_constraint = {}
        
        if 'bond' in constraint:
            bond = constraint['bond']
            new_bond = dict(bond)
            # atom1 and atom2 are FlowList [chain_id, res_idx, atom_name]
            for atom_key in ('atom1', 'atom2'):
                atom = list(bond[atom_key])
                if len(atom) >= 2 and isinstance(atom[1], int):
                    atom[1] = adjust_constraint_residue_index(atom[0], atom[1], offsets)
                new_bond[atom_key] = FlowList(atom)
            new_constraint['bond'] = new_bond
        
        elif 'contact' in constraint:
            contact = constraint['contact']
            new_contact = dict(contact)
            for token_key in ('token1', 'token2'):
                token = list(contact[token_key])
                if len(token) >= 2 and isinstance(token[1], int):
                    token[1] = adjust_constraint_residue_index(token[0], token[1], offsets)
                new_contact[token_key] = FlowList(token)
            new_constraint['contact'] = new_contact
        
        elif 'pocket' in constraint:
            pocket = constraint['pocket']
            new_pocket = dict(pocket)
            # Adjust contacts list: each contact is FlowList [chain_id, res_idx]
            new_contacts = []
            for contact_atom in pocket.get('contacts', []):
                atom = list(contact_atom)
                if len(atom) >= 2 and isinstance(atom[1], int):
                    atom[1] = adjust_constraint_residue_index(atom[0], atom[1], offsets)
                new_contacts.append(FlowList(atom))
            new_pocket['contacts'] = new_contacts
            new_constraint['pocket'] = new_pocket
        
        else:
            new_constraint = constraint
        
        adjusted.append(new_constraint)
    
    return adjusted

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
    compute_affinity = params.get('compute_affinity', True)
    binder_id = params.get('affinity_ligand_id')
    if compute_affinity and binder_id:
        props.append({'affinity': {'binder': binder_id}})
    return props if props else None

def build_constraints(params: Dict) -> Optional[List[Dict]]:
    constraints = []

    # Bond constraints — match pairs of consecutive [chain,res,atom] groups.
    # Multiple bonds are separated by commas. Whitespace between atoms of a bond is tolerated.
    if params.get('constraints_bond'):
        pair_pattern = r'(\[[^\]]+\]\s*\[[^\]]+\])'
        bond_matches = re.findall(pair_pattern, params['constraints_bond'])
        if not bond_matches:
            sys.stderr.write(f"Warning: constraints_bond provided but no valid atom pairs found: '{params['constraints_bond']}'\n")
        for bond_pair in bond_matches:
            atoms = parse_constraint_atoms(bond_pair)
            if len(atoms) == 2:
                constraints.append({
                    'bond': {
                        'atom1': atoms[0],
                        'atom2': atoms[1]
                    }
                })
    
    # Explicit pocket constraints (pocket1 and pocket2)
    for i in range(1, 3):
        binder = params.get(f'constraints_pocket{i}_binder')
        contacts_str = params.get(f'constraints_pocket{i}_contacts')
        if binder and contacts_str:
            # Parse contacts: '[A,53],[A,57],[A,61]' → list of [chain_id, idx] pairs
            contact_pairs = re.findall(r'\[([^\]]+)\]', contacts_str)
            contacts_list = [FlowList(_parse_contact_pair(p)) for p in contact_pairs]
            pocket = {
                'pocket': {
                    'binder': binder,
                    'contacts': contacts_list
                }
            }
            max_dist = params.get(f'constraints_pocket{i}_max_distance')
            if max_dist is not None:
                pocket['pocket']['max_distance'] = float(max_dist)
            force = params.get(f'constraints_pocket{i}_force')
            if force is not None:
                pocket['pocket']['force'] = str(force).lower() == 'true'
            constraints.append(pocket)

    # Extra pocket constraints (string format)
    if params.get('constraints_pocket_extra'):
        constraints.extend(parse_pocket_constraint(params['constraints_pocket_extra']))

    # Contact constraints
    if params.get('constraints_contacts'):
        constraints.extend(parse_contact_constraint(params['constraints_contacts']))

    return constraints if constraints else None

def build_yaml_structure(params: Dict) -> Dict:
    yaml_data = {'sequences': []}
    msa_idx = 0
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
        # Cyclic sequences (Boltz-2 expects 'cyclic: true' inside the sequence type entry)
        cyclic = params.get(f'{seq_key}_cyclic', False)
        if cyclic:
            seq_entry[seq_type]['cyclic'] = True
        # Get start_index for this sequence to adjust modifications
        start_index = int(params.get(f'{seq_key}_start_index', 1) or 1)
        mods = parse_modifications(params.get(f'{seq_key}_mod'))
        if mods:
            if start_index != 1:
                mods = adjust_modification_positions(mods, start_index)
            seq_entry[seq_type]['modifications'] = mods
        # MSA only for protein
        if seq_type == "protein":
            # Check if pre-existing MSA CSV path is provided
            msa_csv_path = params.get(f'{seq_key}_msa_csv')
            if msa_csv_path:
                # Use the provided MSA file path
                seq_entry[seq_type]['msa'] = msa_csv_path
            elif params.get(f'{seq_key}_genmsa') in [None, True]:
                # Generate MSA - reference expected output filename
                pred_name = params.get('pred_name')
                seq_entry[seq_type]['msa'] = f"{pred_name}_msa_{msa_idx}.csv"
                msa_idx += 1
            else:
                seq_entry[seq_type]['msa'] = "empty"
        yaml_data["sequences"].append(seq_entry)

    # Extra sequences
    seqs_extra = parse_sequences_extra(params.get('sequences_extra'))
    yaml_data["sequences"].extend(seqs_extra)

    # Collect residue offsets from explicit params and extra sequences
    # Note: collect_residue_offsets pops _start_index from seqs_extra entries
    offsets = collect_residue_offsets(params, seqs_extra)

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
    constraints = build_constraints(params)
    if constraints:
        # Adjust constraint residue indices if there are offsets
        if offsets:
            constraints = adjust_constraints(constraints, offsets)
        yaml_data['constraints'] = constraints
    templates = parse_templates_from_params(params)
    if templates:
        yaml_data['templates'] = templates
    if properties := build_properties(params):
        yaml_data['properties'] = properties

    # Store residue offsets in the YAML (Boltz-2 ignores unknown top-level fields)
    # This is used by downstream tools (generate_metadata.py, renumber_structures.py)
    if offsets:
        yaml_data['residue_offsets'] = offsets

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
        # Only do this if compute_affinity is enabled (default: true)
        compute_affinity = cleaned_params.get("compute_affinity", True)
        if compute_affinity and not cleaned_params.get("affinity_ligand_id"):
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
        # Only validate affinity if compute_affinity is enabled
        if cleaned_params.get("compute_affinity", True):
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