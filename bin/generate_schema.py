#!/usr/bin/env python3
"""
Generate complete Nextflow schema with 12 sequences and 12 ligands
for Seqera Platform integration - Each sequence and ligand as separate groups
"""

import json
import copy
from typing import Dict, Any

def load_base_schema(filename: str) -> Dict[str, Any]:
    """Load the base schema template"""
    with open(filename, 'r') as f:
        return json.load(f)

def hide_all_properties(properties: Dict[str, Any]):
    for prop in properties.values():
        prop["hidden"] = True
    return properties

def generate_sequence_properties(num: int, is_first: bool = False) -> Dict[str, Any]:
    """Generate properties for a sequence (seq1, seq2, etc.)"""
    
    # Base properties for all sequences
    properties = {
        f"seq{num}_name": {
            "type": "string",
            "description": f"Sequence {num} name/description (optional)",
            "help_text": "Optional name used in plots. Only letters, numbers, underscores, and dashes allowed.",
            "pattern": "^[a-zA-Z0-9_-]+$"
       },
        f"seq{num}_type": {
            "type": "string",
            "description": f"Sequence {num} type",
            "enum": ["protein", "dna", "rna"],
            "help_text": f"Choose the biomolecule type for sequence {num}. Examples: protein, dna, rna"
        },
        f"seq{num}_seq": {
            "type": "string",
            "description": f"Sequence {num} amino acid/nucleotide sequence", 
            "pattern": "^[ACDEFGHIKLMNPQRSTUVWYXZ*-]+$",
            "minLength": 1,
            "maxLength": 10000,
            "help_text": "Single-letter amino acid codes for proteins, ATCG for DNA, AUCG for RNA. Example: ADNKFNKEQQNAFYEILHLPNLTEEQRNGFIQSLKDDPSVSKEILAEAKKLNDAQAPK"
        },
        f"seq{num}_chain_id": {
            "type": "string", 
            "description": f"Unique sequence {num} chain ID(s)",
            "pattern": "^[A-Za-z0-9]+(,[A-Za-z0-9]+)*$",
            "help_text": "Unique single chain ID (e.g., 'A') or multiple IDs separated by commas for multiple copies (e.g., 'A,B,C'). Examples: A, B,C, H,L"
        },
        f"seq{num}_mod": {
            "type": "string",
            "description": f"Sequence {num} modifications e.g. non-canonical amino acids",
            "pattern": "^(\\[\\d+,[A-Z0-9]{3}\\])+$",
            "help_text": "Format: [position,CCD_code][position,CCD_code]... e.g. '[33,SEP]', '[10,PTR][25,TPO]', '[1,ACE][50,NH2]'"
        },
        f"seq{num}_cyclic": {
            "type": "boolean",
            "hidden": True,
            "description": f"Is sequence {num} cyclic?",
            "help_text": "Set to true for cyclic peptides or circular DNA/RNA"
        },
        f"seq{num}_genmsa": {
            "type": "boolean",
            "hidden": True, 
            "description": f"Generate MSA for sequence {num}",
            "default": True,
            "help_text": "Generate multiple sequence alignment to improve prediction accuracy"
        }
    }
    
    # Hide individual parameters for sequences 3-12
    if num > 2:
        return hide_all_properties(properties)
            
    return properties

def generate_ligand_properties(num: int) -> Dict[str, Any]:
    """Generate properties for a ligand (lig1, lig2, etc.)"""
    
    properties = {
        f"lig{num}_type": {
            "type": "string",
            "description": f"Ligand {num} format",
            "enum": ["CCD", "SMILES"],
            "help_text": "CCD: Chemical Component Dictionary code. SMILES: Simplified molecular-input line-entry system. Examples: CCD, SMILES"
        },
        f"lig{num}_chem": {
            "type": "string", 
            "description": f"Ligand {num} chemical identifier",
            "help_text": "CCD code (e.g., 'ATP', 'ZN') or SMILES string (e.g., 'CCO', '[Ca+2]'). Examples: ZN, ATP, CCO, [Ca+2], CC(=O)Nc1ccc(O)cc1"
        },
        f"lig{num}_chain_id": {
            "type": "string",
            "description": f"Unique ligand {num} chain ID(s)",
            "pattern": "^[A-Za-z0-9]+(,[A-Za-z0-9]+)*$", 
            "default": None,
            "help_text": "Unique chain identifier(s) for the ligand e.g. 'B'. Use multiple IDs separated by commas for multiple copies e.g. 'D,E'"
        },
    }
    
    # Hide individual parameters for ligands 3-12
    if num > 2:
        return hide_all_properties(properties)
            
    return properties

def generate_full_schema(base_schema: Dict[str, Any], num_sequences: int = 12, num_ligands: int = 12) -> Dict[str, Any]:
    schema = copy.deepcopy(base_schema)
    new_defs = {}
    allof_refs = []
    # Always add IO and prediction options first
    for sec in ["input_output_options", "prediction_options"]:
        new_defs[sec] = schema["definitions"][sec]
        allof_refs.append({"$ref": f"#/definitions/{sec}"})
    # Sequence and ligand groups
    for i in range(1, num_sequences + 1):
        new_defs[f"sequence_{i}_options"] = {
            "title": f"Sequence {i} options",
            "type": "object",
            "fa_icon": "fas fa-dna",
            "description": f"Configure sequence {i} (protein, DNA, or RNA)",
            "properties": generate_sequence_properties(i),
            **({"required": [f"seq{i}_type", f"seq{i}_chain_id", f"seq{i}_seq"]} if i == 1 else {}),
            **({"hidden": True} if i > 2 else {})
        }
        allof_refs.append({"$ref": f"#/definitions/sequence_{i}_options"})
    for i in range(1, num_ligands + 1):
        new_defs[f"ligand_{i}_options"] = {
            "title": f"Ligand {i} options",
            "type": "object",
            "fa_icon": "fas fa-prescription-bottle",
            "description": f"Configure ligand {i}",
            "properties": generate_ligand_properties(i),
            **({"hidden": True} if i > 2 else {})
        }
        allof_refs.append({"$ref": f"#/definitions/ligand_{i}_options"})
    # Other sections (constraint and compute)
    for sec in ["affinity_options","constraint_options","template_options","compute_options"]:
        if sec in schema["definitions"]:
            new_defs[sec] = schema["definitions"][sec]
            allof_refs.append({"$ref": f"#/definitions/{sec}"})
    schema["definitions"] = new_defs
    schema["allOf"] = allof_refs
    return schema

def main():
    """Main function"""
    print("Loading base schema...")
    base_schema = load_base_schema("nextflow_base_schema.json")
    
    print("Generating full schema with 12 individual sequence and ligand groups...")
    full_schema = generate_full_schema(base_schema)
    
    print("Writing output schema...")
    with open("nextflow_schema.json", "w") as f:
        json.dump(full_schema, f, indent=2)
    
    print("✅ Schema generated successfully!")
    print("📋 Configuration:")
    print("   - 12 individual sequence groups (Sequence 1-12 options)")
    print("   - 12 individual ligand groups (Ligand 1-12 options)")
    print("   - seq1 is required with defaults")
    print("   - Sequences 3-12 groups AND individual parameters are hidden")
    print("   - Ligands 3-12 groups AND individual parameters are hidden")
    print("   - Constraint options are hidden")
    print("   - Only boltz_diffusion_samples is visible in Boltz options")
    print("   - All compute options are hidden")

if __name__ == "__main__":
    main()
