#!/usr/bin/env python3
"""
Generate Nextflow schema for nf-boltz with 4 sequences, 4 ligands,
plus 'extra' blocks for sequences and ligands (in string 'block' format).
"""

import json
import copy
from typing import Dict, Any

def load_base_schema(filename: str) -> Dict[str, Any]:
    with open(filename, 'r') as f:
        return json.load(f)

def generate_sequence_properties(num: int) -> Dict[str, Any]:
    return {
        f"seq{num}_type": {
            "type": "string",
            "enum": ["protein", "dna", "rna"],
            "description": f"Sequence {num}: type",
            "help_text": "Choose one: protein, dna, or rna."
        },
        f"seq{num}_seq": {
            "type": "string",
            "description": f"Sequence {num}: amino acid/nucleotide sequence",
            "minLength": 1,
            "maxLength": 5000,
            "help_text": ("Sequence using single-letter amino acid codes for proteins, ATCG for DNA/ RNA. "
                          "Type will be auto-detected: do not prefix with anything. E.g., 'MDYK...'."),
            "pattern": "^[A-Za-z*.,\\-\\s]+$"  # loose: actual type/chars validated downstream
        },
        f"seq{num}_chain_id": {
            "type": "string",
            "description": f"Sequence {num}: unique chain ID(s)",
            "pattern": "^[A-Za-z0-9_]+(,[A-Za-z0-9_]+)*$",
            "help_text": "Unique single chain id ('A') or comma-separated list for multiple copies ('A,B,C')."
        },
        f"seq{num}_mod": {
            "type": "string",
            "hidden": True,
            "description": f"Sequence {num}: modifications (optional)",
            "pattern": "^(\\[\\d+,([A-Z0-9]{3})\\])*$",
            "help_text": "Format: [pos,CCD][pos,CCD]... (e.g., '[10,SEP][55,PTR]'), or leave blank if none."
        },
        f"seq{num}_cyclic": {
            "type": "boolean",
            "hidden": True,
            "description": f"Is sequence {num} cyclic?",
            "help_text": "Set to true for cyclic peptides/DNA/RNA.",
            "default": False
        },
        f"seq{num}_genmsa": {
            "type": "boolean",
            "hidden": True,
            "description": f"Generate MSA for sequence {num} (if protein)?",
            "default": True,
            "help_text": "Generate multiple sequence alignment to improve accuracy (protein only)"
        }
    }

def generate_ligand_properties(num: int) -> Dict[str, Any]:
    return {
        f"lig{num}_chain_id": {
            "type": "string",
            "description": f"Ligand {num}: unique chain ID(s)",
            "pattern": "^[A-Za-z0-9_]+(,[A-Za-z0-9_]+)*$",
            "help_text": "Single id ('D') or comma-separated ('G,H')."
        },
        f"lig{num}_ccd": {
            "type": "string",
            "description": f"Ligand {num}: CCD code",
            "help_text": "E.g. 'ZN', 'ATP'. Leave blank if using only SMILES."
        },
        f"lig{num}_smiles": {
            "type": "string",
            "description": f"Ligand {num}: SMILES string",
            "help_text": "E.g. '[Ca+2]', 'CCO', 'C1CCCCC1'. If provided, CCD is ignored."
        }
    }

def generate_full_schema(base_schema: Dict[str, Any]) -> Dict[str, Any]:
    schema = copy.deepcopy(base_schema)
    new_defs = {}
    allof_refs = []

    # Always add core sections
    for sec in ["input_output_options", "prediction_options"]:
        new_defs[sec] = schema["definitions"][sec]
        allof_refs.append({"$ref": f"#/definitions/{sec}"})

    # 4 explicit sequence blocks
    for i in range(1, 5):
        new_defs[f"sequence_{i}_options"] = {
            "title": f"Sequence {i} options",
            "type": "object",
            "fa_icon": "fas fa-dna",
            "description": f"Configure sequence {i} (protein, DNA or RNA chain)",
            "properties": generate_sequence_properties(i)
        }
        allof_refs.append({"$ref": f"#/definitions/sequence_{i}_options"})
    # --- sequences_extra ---
    new_defs["sequences_extra_options"] = {
        "title": "Extra sequences (advanced)",
        "type": "object",
        "fa_icon": "fas fa-plus-square",
        "description": "Add additional input sequences here. Format: space separated blocks: [ID(s)][SEQUENCE] [ID(s)][SEQUENCE] ...",
        "properties": {
            "sequences_extra": {
                "type": "string",
                "default": None,
                "description": "Extra sequences block-string",
                "pattern": "^(\\[[A-Za-z0-9_,]+\\]\\[[A-Za-z*.,\\-\\s]+\\]\\s*)*$",
                "help_text": (
                    "Extra sequences in block format: [IDs][SEQUENCE] [IDs][SEQUENCE] ... "
                    "IDs = comma separated for multi-chain. Sequence type is auto-detected. Example: '[E][MAAAD] [F,G][ACTGACTG]'. "
                    "No modifications or cyclic/msa control."
                ),
                "hidden": True  # Hide from default UI
            }
        }
    }
    allof_refs.append({"$ref": "#/definitions/sequences_extra_options"})

    # 4 explicit ligand blocks
    for i in range(1, 5):
        new_defs[f"ligand_{i}_options"] = {
            "title": f"Ligand {i} options",
            "type": "object",
            "fa_icon": "fas fa-prescription-bottle",
            "description": f"Configure ligand {i} (small molecule/ion)",
            "properties": generate_ligand_properties(i)
        }
        allof_refs.append({"$ref": f"#/definitions/ligand_{i}_options"})


    # --- ligands_extra ---
    new_defs["ligands_extra_options"] = {
        "title": "Extra ligands (advanced)",
        "type": "object",
        "fa_icon": "fas fa-plus-square",
        "description": "Add additional ligands here. Format: space separated blocks: [ID(s)][CCD][SMILES] [ID(s)][CCD][SMILES] ...",
        "properties": {
            "ligands_extra": {
                "type": "string",
                "default": None,
                "description": "Extra ligands block-string",
                "pattern": (
                    r"^(\[[A-Za-z0-9_,]+\]"                # IDs, comma
                    r"\[[A-Za-z0-9\-+]*\]"                 # CCD (may be blank)
                    r"\[[^\[\]]*\]\s*)*$"                  # SMILES (may be blank, allow anything but brackets inside)
                ),
                "help_text": (
                    "Extra ligands: space separated blocks: [IDs][CCD][SMILES]. "
                    "IDs: comma-separated for multiple copies. For SMILES only, leave CCD blank '[E][][C1CCCCC1]'; for CCD only, leave SMILES blank. "
                    "Example: '[H][ZN][] [I][][N[C@@H](Cc1ccc(O)cc1)C(=O)O]'"
                ),
                "hidden": True  # Hide from UI
            }
        }
    }
    allof_refs.append({"$ref": "#/definitions/ligands_extra_options"})

    # Add other sections (constraints, templates, affinity, compute)
    for sec in ["affinity_options", "constraint_options", "template_options", "compute_options"]:
        if sec in schema["definitions"]:
            new_defs[sec] = schema["definitions"][sec]
            allof_refs.append({"$ref": f"#/definitions/{sec}"})

    schema["definitions"] = new_defs
    schema["allOf"] = allof_refs
    return schema

def main():
    base_schema = load_base_schema("nextflow_base_schema.json")
    full_schema = generate_full_schema(base_schema)
    with open("nextflow_schema.json", "w") as f:
        json.dump(full_schema, f, indent=2)
    print("✅ Nextflow schema with 4 sequences/ligands + extra fields generated successfully!")

if __name__ == "__main__":
    main()
