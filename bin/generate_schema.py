#!/usr/bin/env python3
"""
Generate complete Nextflow schema with 12 sequences and 12 ligands
for Seqera Platform integration
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
            "help_text": f"Choose the biomolecule type for sequence {num}"
        },
        f"seq{num}_id": {
            "type": "string", 
            "description": f"Sequence {num} chain ID(s)",
            "pattern": "^[A-Za-z0-9]+(,[A-Za-z0-9]+)*$",
            "help_text": "Single chain ID (e.g., 'A') or multiple IDs separated by commas (e.g., 'A,B,C')",
            "examples": ["A", "B,C", "H,L"]
        },
        f"seq{num}_seq": {
            "type": "string",
            "description": f"Sequence {num} amino acid/nucleotide sequence", 
            "pattern": "^[ACDEFGHIKLMNPQRSTVWYXZ*-]+$",
            "minLength": 1,
            "maxLength": 10000,
            "help_text": "Single-letter amino acid codes for proteins, ATCG for DNA, AUCG for RNA"
        },
        f"seq{num}_mod": {
            "type": "string",
            "description": f"Sequence {num} modifications",
            "pattern": "^(\\[\\d+,[A-Z]{3}\\])(,\\[\\d+,[A-Z]{3}\\])*$",
            "help_text": "Format: [position,CCD_code]. Multiple modifications separated by commas.",
            "examples": ["[33,SEP]", "[10,PTR],[25,TPO]", "[1,ACE],[50,NH2]"]
        },
        f"seq{num}_cyclic": {
            "type": "boolean",
            "description": f"Is sequence {num} cyclic?",
            "default": False,
            "help_text": "Set to true for cyclic peptides or circular DNA/RNA"
        },
        f"seq{num}_genmsa": {
            "type": "boolean", 
            "description": f"Generate MSA for sequence {num}",
            "default": False,
            "help_text": "Generate multiple sequence alignment to improve prediction accuracy"
        }
    }
    
    # First sequence specific settings
    if is_first:
        properties[f"seq{num}_type"]["default"] = "protein"
        properties[f"seq{num}_id"]["default"] = "A"
        properties[f"seq{num}_seq"]["default"] = "ADNKFNKEQQNAFYEILHLPNLTEEQRNGFIQSLKDDPSVSKEILAEAKKLNDAQAPK"
        properties[f"seq{num}_mod"]["default"] = "[33,SEP],[42,SEP]"
        properties[f"seq{num}_genmsa"]["default"] = True
    
    # Hide sequences 4-12 by default
    if num > 3:
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
            "default": None,  # All ligands default to null
            "help_text": "Chain identifier(s) for the ligand",
            "examples": ["D", "L1,L2"]
        },
        f"lig{num}_type": {
            "type": "string",
            "description": f"Ligand {num} format",
            "enum": ["CCD", "SMILES"],
            "default": None,  # All ligands default to null
            "help_text": "CCD: Chemical Component Dictionary code. SMILES: Simplified molecular-input line-entry system"
        },
        f"lig{num}_chem": {
            "type": "string", 
            "description": f"Ligand {num} chemical identifier",
            "default": None,  # All ligands default to null
            "help_text": "CCD code (e.g., 'ATP', 'ZN') or SMILES string (e.g., 'CCO', '[Ca+2]')",
            "examples": ["ZN", "ATP", "CCO", "[Ca+2]", "CC(=O)Nc1ccc(O)cc1"]
        }
    }
    
    # Hide ligands 4-12 by default
    if num > 3:
        for key in properties:
            properties[key]["hidden"] = True
            
    return properties

def generate_full_schema(base_schema: Dict[str, Any]) -> Dict[str, Any]:
    """Generate the complete schema with all sequences and ligands"""
    
    schema = copy.deepcopy(base_schema)
    
    # Create sequence options section
    seq_properties = {}
    for i in range(1, 13):
        seq_props = generate_sequence_properties(i, is_first=(i==1))
        seq_properties.update(seq_props)
    
    schema["definitions"]["sequence_options"] = {
        "title": "Sequence options",
        "type": "object",
        "fa_icon": "fas fa-dna", 
        "description": "Configure protein, DNA, and RNA sequences",
        "required": ["seq1_type", "seq1_id", "seq1_seq"],
        "properties": seq_properties
    }
    
    # Create ligand options section
    lig_properties = {}
    for i in range(1, 13):
        lig_props = generate_ligand_properties(i)
        lig_properties.update(lig_props)
        
    schema["definitions"]["ligand_options"] = {
        "title": "Ligand options",
        "type": "object", 
        "fa_icon": "fas fa-prescription-bottle",
        "description": "Configure small molecules, ions, and other ligands",
        "properties": lig_properties
    }
    
    # Add constraint defaults
    schema["definitions"]["constraint_options"]["properties"]["constraints_bond"]["default"] = "[A,33,O][B,1,N] [B,1,CA][C,5,CA]"
    schema["definitions"]["constraint_options"]["properties"]["pocket_binder_id"]["default"] = "E" 
    schema["definitions"]["constraint_options"]["properties"]["pocket_contacts"]["default"] = "[A,33][A,34]"
    schema["definitions"]["constraint_options"]["hidden"] = True
    
    # Set up allOf references in the desired order
    schema["allOf"] = [
        {"$ref": "#/definitions/input_output_options"},    # 1. Input/Output options
        {"$ref": "#/definitions/sequence_options"},        # 2. Sequence options  
        {"$ref": "#/definitions/ligand_options"},          # 3. Ligand options
        {"$ref": "#/definitions/constraint_options"},      # 4. Constraint options
        {"$ref": "#/definitions/boltz_options"},           # 5. Boltz-1 prediction options
        {"$ref": "#/definitions/compute_options"}          # 6. Compute options
    ]
    
    # Clean up templates
    templates_to_remove = [k for k in schema["definitions"].keys() if k.endswith("_template")]
    for template in templates_to_remove:
        del schema["definitions"][template]
    
    return schema

def main():
    """Main function"""
    print("Loading base schema...")
    base_schema = load_base_schema("nextflow_base_schema.json")
    
    print("Generating full schema with 12 sequences and 12 ligands...")
    full_schema = generate_full_schema(base_schema)
    
    print("Writing output schema...")
    with open("nextflow_schema.json", "w") as f:
        json.dump(full_schema, f, indent=2)
    
    print("✅ Schema generated successfully!")
    print("📋 Configuration:")
    print("   - 12 sequences (seq1-seq12)")
    print("   - 12 ligands (lig1-lig12) - ALL DEFAULT TO NULL")
    print("   - seq1 is required with defaults")
    print("   - seq4-seq12 and lig4-lig12 are hidden")
    print("   - Only boltz_diffusion_samples is visible in Boltz options")
    print("   - All compute options are hidden")

if __name__ == "__main__":
    main()
