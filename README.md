# nf-boltz

A robust, scalable pipeline to run **Boltz-2** or **Boltz-2** structure predictions with advanced complex modeling, ligand support, and structural constraints, all configured using [Nextflow](https://www.nextflow.io/).

## Introduction

**nf-boltz** orchestrates Boltz structure prediction, handling input preparation, resource management, sequences/ligands complexity, and advanced configuration for constraints and templates.

- Fast setup of complex protein, DNA, and RNA systems.
- Simple parameter specification; flexible for power users via "extra" fields.
- Supports ligands (by CCD or SMILES), pockets, contacts, templates, and batch runs.

Upstream Boltz: https://github.com/jwohlwend/boltz

## Quick Start

Install Nextflow as per their docs.

Clone this repository and run:

    nextflow run main.nf --pred_name my_example --outdir ./results

> By default, the system expects sequences/ligands defined via pipeline parameters.
> Alternatively, you may provide a completed YAML file via `--input_yaml` (see [YAML input mode](#yaml-input-mode) below).

## Input Parameters

All parameters can be provided as Nextflow CLI flags (e.g. -params-file, or individually with --X value).

### Main Input/Output

| Parameter   | Description                                            | Default              |
|-------------|--------------------------------------------------------|----------------------|
| pred_name   | Run name (letters, numbers, _, -)                      | test                 |
| outdir      | Output directory                                       | /vast/scratch/users/hardy.j/boltztest |
| input_yaml  | Optional: YAML file/dir with pre-configured input      | null                 |

### Sequences (Proteins / DNA / RNA)

Supports **up to 4 main ("explicit") sequences** plus additional ones using "extra" string block format.

For each of 4 explicit entries (`seq1`...`seq4`):

- `seqX_type`      — Must be `protein`, `dna`, or `rna`
- `seqX_chain_id`  — Chain ID(s), e.g. 'A' or 'A,B'
- `seqX_seq`       — Single-letter sequence
- `seqX_mod`       — Optional, modifications as `[RESIDX,CCD]...` blocks, e.g. `[10,SEP][55,PTR]`
- `seqX_cyclic`    — true if cyclic (rare; usually false)
- `seqX_genmsa`    — true for MSA (default for protein)

Example:
    seq1_type      = 'protein'
    seq1_chain_id  = 'A'
    seq1_seq       = 'MADSEQYKLV...'

    seq2_type      = 'dna'
    seq2_chain_id  = 'B'
    seq2_seq       = 'ACTGGTTCGA'

For 5 or more sequences, use sequences_extra:
    sequences_extra = '[E][MAAADTESTSEQ] [F,G][CCGUUUAA]'

- Format: [ID(s)][SEQUENCE]
- Chain IDs can be comma-separated for multiple chains.
- Sequence type is auto-detected for extras.

### Ligands (Small Molecules / Ions)

Supports **up to 4 main ligands** plus more using "extra" string block format.

Each of 4 explicit ligands (`lig1`...`lig4`):

- `ligX_chain_id` — Ligand chain IDs (e.g. 'D' or 'D,E')
- `ligX_ccd`      — CCD code (e.g. 'ZN', 'ATP') (optional)
- `ligX_smiles`   — SMILES string (optional; takes precedence over ccd if supplied)

Example:
    lig1_chain_id = 'C,D'
    lig1_ccd      = 'SAH'
    lig1_smiles   = 'CCO'

For 5 or more ligands, use ligands_extra:
    ligands_extra = '[E][ZN][] [F,G][][C1CCCCC1]'

- Format: [ID(s)][CCD][SMILES]
- Either CCD or SMILES required; chain IDs can be comma-separated.

### Constraints (Optional)

#### Bond Constraints

- Format: [CHAIN_ID,RES,ATOM][CHAIN_ID,RES,ATOM] ...
- Example: constraints_bond = '[A,33,O][B,1,N] [H,95,SG][L,23,SG]'

#### Contact Constraints

- Format: Each contact as [CHAIN_ID,IDX/ATOM][CHAIN_ID,IDX/ATOM][DIST][FORCE]
- Example: constraints_contacts = '[A,53][A,57] [A,101][A,130][true]'

#### Pocket Constraints

- Two explicit pockets: constraints_pocket1_binder, constraints_pocket1_contacts, etc.
- More than 2: Use constraints_pocket_extra as space-separated blocks:
  constraints_pocket_extra = "[D][[A,5],[B,12]] [E][[B,99],[C,1]][true]"

### Templates (Optional)

- Up to 2 explicit templates (template1_path, etc.)
- More using template_extra as:
  [PATH][PREDICTION_IDS][TEMPLATE_IDS][THRESHOLD] [PATH2][...]

### Affinity (Optional)

- affinity_ligand_id — Ligand chain ID to use for affinity calculations.
- Additional control over affinity prediction via advanced parameters.

### Compute and Advanced

- Resources (boltz2_gpu_model, boltz2_cpus, memory, etc.)
- Output format (output_format = 'mmcif' by default).
- Extra Boltz args (extra_config) if needed.

## YAML Input Mode

You can supply a complete, prepared input YAML using --input_yaml.
> In this case, *all other sequence/ligand/constraint/template parameters are ignored*.

For YAML format details, see Boltz documentation:
https://github.com/jwohlwend/boltz/blob/main/docs/prediction.md

## Examples

#### Minimal Protein Example

    nextflow run main.nf \
      --pred_name myprot \
      --seq1_type protein \
      --seq1_chain_id A \
      --seq1_seq MVLSPADKTNVKAA...

#### Protein-Ligand Complex (With SMILES)

    nextflow run main.nf \
      --pred_name myligand \
      --seq1_type protein \
      --seq1_chain_id A \
      --seq1_seq MSEQN... \
      --lig1_chain_id E \
      --lig1_smiles '[Na+1]'

#### Complex with Extras

    nextflow run main.nf \
      --pred_name batch \
      --seq1_type protein --seq1_chain_id A --seq1_seq MA... \
      --seq2_type dna     --seq2_chain_id B --seq2_seq ACTGA... \
      --sequences_extra '[F][MGKQHVPTAE...] [K,L][TTAGGG]' \
      --ligands_extra '[X,Y][HEM][] [Z][][C1CCCCC1]'

## Output

- Structure predictions in output directory (PDB/mmCIF).
- JSON or YAML summary of job configuration.
- Optional per-model confidence/PAE/PDE as selected.

## Tips & Troubleshooting

- Explicit sequence types: For seq1–seq4, always specify seqX_type. "Extras" will be auto-detected.
- Chain IDs: must be unique across all sequences and ligands.
- Ligands: Provide SMILES if possible. If both CCD and SMILES are given, SMILES is used.
- Constraints and advanced features: See Boltz constraints documentation and pipeline comments.
- For large/complex jobs: Use -with-dag or -with-report for Nextflow run visualization/reporting.

## References

- [Main Boltz repository](https://github.com/jwohlwend/boltz) [bioRxiv](https://doi.org/10.1101/2025.06.14.659707)
- [Low-memory Boltz-2 implementation](https://github.com/tlitfin/lmi4boltz) [bioRxiv](https://doi.org/10.1101/2025.10.29.684571)
- [Boltz prediction input format](https://github.com/jwohlwend/boltz/blob/main/docs/prediction.md)
- [Nextflow documentation](https://www.nextflow.io/)

## Credits

**nf-boltz** was originally written by [Joshua Hardy](mailto:hardy.j@wehi.edu.au).
