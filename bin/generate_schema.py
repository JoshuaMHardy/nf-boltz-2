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

def generate_sequence_properties(num: int, is_first: bool = False) -> Dict[str, Any]:
    """Generate properties for a sequence (seq1, seq2, etc.)"""
    
    # Base properties for all sequences
    properties = {
        f"seq{num}_type": {
            "type": "string",
            "description": f"Sequence {num} type",
            "enum": ["protein", "dna", "rna"],
            "help_text": f"Choose the biomolecule type for sequence {num}. Examples: protein, dna, rna"
        },
        f"seq{num}_id": {
            "type": "string", 
            "description": f"Sequence {num} chain ID(s)",
            "pattern": "^[A-Za-z0-9]+(,[A-Za-z0-9]+)*$",
            "help_text": "Single chain ID (e.g., 'A') or multiple IDs separated by commas (e.g., 'A,B,C'). Examples: A, B,C, H,L"
        },
        f"seq{num}_seq": {
            "type": "string",
            "description": f"Sequence {num} amino acid/nucleotide sequence", 
            "pattern": "^[ACDEFGHIKLMNPQRSTVWYXZ*-]+$",
            "minLength": 1,
            "maxLength": 10000,
            "help_text": "Single-letter amino acid codes for proteins, ATCG for DNA, AUCG for RNA. Example: ADNKFNKEQQNAFYEILHLPNLTEEQRNGFIQSLKDDPSVSKEILAEAKKLNDAQAPK"
        },
        f"seq{num}_mod": {
            "type": "string",
            "description": f"Sequence {num} modifications",
            "pattern": "^(\\[\\d+,[A-Z0-9]{3}\\])+$",
            "help_text": "Format: [position,CCD_code][position,CCD_code]... e.g. '[33,SEP]', '[10,PTR][25,TPO]', '[1,ACE][50,NH2]'"
        },
        f"seq{num}_cyclic": {
            "type": "boolean",
            "description": f"Is sequence {num} cyclic?",
            "help_text": "Set to true for cyclic peptides or circular DNA/RNA"
        },
        f"seq{num}_genmsa": {
            "type": "boolean", 
            "description": f"Generate MSA for sequence {num}",
            "default": True,
            "help_text": "Generate multiple sequence alignment to improve prediction accuracy"
        }
    }
    
    # First sequence specific settings
    if is_first:
        properties[f"seq{num}_type"]["default"] = "protein"
        properties[f"seq{num}_id"]["default"] = "A"
        properties[f"seq{num}_seq"]["default"] = "ADNKFNKEQQNAFYEILHLPNLTEEQRNGFIQSLKDDPSVSKEILAEAKKLNDAQAPK"
        properties[f"seq{num}_mod"]["default"] = "[33,SEP][42,SEP]"
        properties[f"seq{num}_genmsa"]["default"] = True
    
    # Hide individual parameters for sequences 3-12
    if num > 2:
        for key in properties:
            properties[key]["hidden"] = True
            
    return properties

def generate_ligand_properties(num: int) -> Dict[str, Any]:
    """Generate properties for a ligand (lig1, lig2, etc.)"""
    
    properties = {
        f"lig{num}_id": {
            "type": "string",
            "description": f"Ligand {num} chain ID(s)",
            "pattern": "^[A-Za-z0-9]+(,[A-Za-z0-9]+)*$", 
            "default": None,
            "help_text": "Chain identifier(s) for the ligand. Examples: D, L1,L2"
        },
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
        }
    }
    
    # Hide individual parameters for ligands 3-12
    if num > 2:
        for key in properties:
            properties[key]["hidden"] = True
            
    return properties

def generate_full_schema(base_schema: Dict[str, Any]) -> Dict[str, Any]:
    """Generate the complete schema with all sequences and ligands"""
    
    schema = copy.deepcopy(base_schema)
    
    # Start with base definitions
    new_definitions = {
        "input_output_options": schema["definitions"]["input_output_options"],
        "boltz_options": schema["definitions"]["boltz_options"]
    }
    
    # Create individual sequence groups
    sequence_refs = []
    for i in range(1, 13):
        seq_props = generate_sequence_properties(i, is_first=(i==1))
        
        group_definition = {
            "title": f"Sequence {i} options",
            "type": "object",
            "fa_icon": "fas fa-dna",
            "description": f"Configure sequence {i} (protein, DNA, or RNA)",
            "properties": seq_props
        }
        
        # Add required fields for sequence 1
        if i == 1:
            group_definition["required"] = ["seq1_type", "seq1_id", "seq1_seq"]
        
        # Hide sequences 3-12
        if i > 2:
            group_definition["hidden"] = True
            
        new_definitions[f"sequence_{i}_options"] = group_definition
        sequence_refs.append({"$ref": f"#/definitions/sequence_{i}_options"})
    
    # Create individual ligand groups
    ligand_refs = []
    for i in range(1, 13):
        lig_props = generate_ligand_properties(i)
        
        group_definition = {
            "title": f"Ligand {i} options",
            "type": "object",
            "fa_icon": "fas fa-prescription-bottle",
            "description": f"Configure ligand {i} (small molecules, ions, etc.)",
            "properties": lig_props
        }
        
        # Hide ligands 3-12
        if i > 2:
            group_definition["hidden"] = True
            
        new_definitions[f"ligand_{i}_options"] = group_definition
        ligand_refs.append({"$ref": f"#/definitions/ligand_{i}_options"})
    
    # Add constraint and compute options
    constraint_options = copy.deepcopy(schema["definitions"]["constraint_options"])
    constraint_options["hidden"] = True
    
    new_definitions["constraint_options"] = constraint_options
    new_definitions["compute_options"] = schema["definitions"]["compute_options"]
    
    # Update schema with new definitions
    schema["definitions"] = new_definitions
    
    # Set up allOf references in the desired order
    allof_refs = [
        {"$ref": "#/definitions/input_output_options"},
        {"$ref": "#/definitions/boltz_options"}
    ]
    
    # Add sequence references
    allof_refs.extend(sequence_refs)
    
    # Add ligand references  
    allof_refs.extend(ligand_refs)
    
    # Add constraint and compute
    allof_refs.extend([
        {"$ref": "#/definitions/constraint_options"},
        {"$ref": "#/definitions/compute_options"}
    ])
    
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
