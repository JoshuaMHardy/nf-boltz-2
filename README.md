# nf-boltz-2

A robust, scalable pipeline to run **Boltz-2** structure predictions with advanced complex modeling, ligand support, and structural constraints, with local multiple sequence alignment (MSA) using [Nextflow](https://www.nextflow.io/).

<img width="800" src="img/nf-boltz-2-workflow.png">

## Introduction

**nf-boltz-2** orchestrates Boltz structure prediction, handling input preparation, resource management, sequences/ligands complexity, and advanced configuration for constraints and templates.

- Fast setup of complex protein, DNA, and RNA systems.
- Simple parameter specification; flexible for power users via "extra" fields.
- Supports ligands (by CCD or SMILES), pockets, contacts, templates, and batch runs.

If you find this useful in your research please cite this [workflow](https://doi.org/10.48546/workflowhub.workflow.2318.1) and the open-source tools that make workflows like this possible (see [Citations](#Citations)).

## Setup

This pipeline runs every step in an Apptainer/Singularity container - nothing needs to be
installed on the host besides [Nextflow](https://www.nextflow.io/) and Apptainer itself.

### 1. Build the containers

```
apptainer build boltz2.sif apptainer/boltz2.def
apptainer build boltz2-lowmem.sif apptainer/boltz2-lowmem.def   # optional, for large/low-memory GPU jobs
apptainer build colabfold.sif apptainer/colabfold.def
apptainer build python_tools.sif apptainer/python_tools.def
```

ChimeraX is free for non-commercial/academic use but its `.deb` cannot be redistributed here.
Agree to the license and [download ChimeraX for Ubuntu 24.04](https://www.cgl.ucsf.edu/chimerax/download.html)
into `apptainer/`, then build:
```
apptainer build chimerax.sif apptainer/chimerax.def
```
Set `--enable_visualisation false` if you don't want to build/use this container.

By default the pipeline looks for these `.sif` files in the repository root
(`boltz2_container`, `colabfold_container`, etc. in [nextflow.config](nextflow.config)).
Override with `--boltz2_container /path/to/boltz2.sif` etc. if you keep them elsewhere.

### 2. Boltz-2 model weights

`boltz predict --cache <dir>` downloads model weights automatically into `boltz2_models` on
first run if the directory is empty, so no manual step is required - just make sure the path
(default `./boltzmodels`) is writable and has network access on first run.

### 3. ColabFold sequence databases (for MSA generation)

MSA generation runs `colabfold_search` against local MMseqs2 databases. Edit
[setup_colabfold_databases.sbatch](setup_colabfold_databases.sbatch) (`WORK_DIR`, `CONTAINER`,
and SLURM resource directives) for your cluster, then submit it. The databases are large
(several hundred GB) and download/indexing can take many hours. Point `--colabfold_dbs` at
the resulting directory (default `./colabfold_dbs`).

If you don't need MSAs (e.g. `msa: empty` inputs only, or all `seqX_genmsa = false`), you can
skip this step entirely.

### 4. Execution profile

The pipeline ships without an institution-specific executor/scheduler config by default -
pick a profile for your environment:

| Profile      | Use case                                                                 |
|--------------|---------------------------------------------------------------------------|
| `standalone` | Single machine / workstation with a local GPU (local executor, enables Apptainer) |
| `slurm`      | Generic SLURM cluster (enables Apptainer + SLURM executor); combine with `--gpu_queue`, `--boltz2_container`, etc. for your site |
| `milton`     | WEHI's Milton SLURM cluster (institutional container/database paths + GPU queue name) |

For other SLURM clusters, use `-profile slurm` and either pass container/database paths as
CLI flags, or copy the `milton` profile block in [nextflow.config](nextflow.config) as a
template for your own site-specific `params` overrides.

## Quick Start

Clone this repository, complete [Setup](#setup) above, then run:

    nextflow run main.nf --pred_name my_example --outdir ./results -profile standalone

> By default, the system expects sequences/ligands defined via pipeline parameters.
> Alternatively, you may provide a completed YAML file via `--input_yaml` (see [YAML input mode](#yaml-input-mode) below).

For more information on Nextflow setup and use see the official [Nextflow documentation](https://www.nextflow.io/).

## Input Parameters

All parameters can be provided as Nextflow CLI flags (e.g. -params-file, or individually with --X value).

### Main Input/Output

| Parameter         | Description                                                                                                                | Default |
|-------------------|-----------------------------------------------------------------------------------------------------------------------------|---------|
| pred_name         | Run name (letters, numbers, _, -)                                                                                          | test    |
| outdir            | Output directory (relative or absolute path)                                                                               | boltz2  |
| use_outdir_prefix | Whether to prepend `outdir_prefix` to `outdir`. Set `true` when `outdir_prefix` is used, `false` for a plain relative/absolute `outdir` | false   |
| outdir_prefix     | Optional output directory prefix, e.g. `/vast/scratch/users/$USER/` (supports automatic username substitution on Seqera)   | null    |
| input_yaml        | Optional: YAML file/dir with pre-configured input. If set, all other sequence/ligand/constraint/template parameters below are ignored | null    |

### Sequences (Proteins / DNA / RNA)

Supports **up to 4 main ("explicit") sequences** plus additional ones using "extra" string block format.

For each of 4 explicit entries (`seq1`...`seq4`):

- `seqX_type`       — Must be `protein`, `dna`, or `rna`
- `seqX_chain_id`   — Chain ID(s), e.g. 'A' or 'A,B'
- `seqX_seq`        — Single-letter sequence
- `seqX_mod`        — Optional, modifications as `[RESIDX,CCD]` blocks, comma-separated, e.g. `[10,SEP],[55,PTR]`
- `seqX_cyclic`     — true if cyclic (rare; usually false). Emits `cyclic: true` in the Boltz-2 YAML.
- `seqX_genmsa`     — true for MSA generation (default for protein); false skips MSA (`msa: empty`)
- `seqX_msa_csv`    — Path to a pre-existing MSA CSV file, skips MSA generation for this sequence
- `seqX_start_index`— Starting residue index (default: 1) for partial sequences; modifications and constraints use real biological numbering and are auto-adjusted for Boltz-2, then converted back in outputs/plots

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
- SMILES containing brackets are supported (bracket-depth parsed, not regex), e.g.
  `'[D][][[Na+]],[E][][CC(=O)[O-]]'`.

### Constraints (Optional)

#### Bond Constraints

- Format: [CHAIN_ID,RES,ATOM][CHAIN_ID,RES,ATOM] ...
- Example: constraints_bond = '[A,33,O][B,1,N] [H,95,SG][L,23,SG]'
- Optional whitespace between the two atoms of a pair is tolerated.

#### Contact Constraints

- Format: Each contact as [CHAIN_ID,IDX/ATOM][CHAIN_ID,IDX/ATOM][DIST][FORCE]
- Example: constraints_contacts = '[A,53][A,57] [A,101][A,130][true]'

#### Pocket Constraints

- Two explicit pockets: constraints_pocket1_binder, constraints_pocket1_contacts, etc.
- More than 2: Use constraints_pocket_extra as space-separated blocks:
  constraints_pocket_extra = "[D][[A,5],[B,12]] [E][[B,99],[C,1]][true]"

### Templates (Optional)

Supports **up to 2 explicit templates** (`template1`...`template2`) plus additional ones using an "extra" string block format.

For each of 2 explicit entries:

- `templateX_pdb_code`      — 4-character PDB code; downloads the CIF file from RCSB
- `templateX_path`         — Path to a local CIF file (alternative to `pdb_code`)
- `templateX_template_id`  — Template chain ID(s), comma-separated. Leave blank to auto-match
- `templateX_prediction_id`— Prediction chain ID(s), comma-separated. Leave blank to auto-match
- `templateX_threshold`    — Optional max RMSD (Å) between template and prediction, enforced with a force potential. If omitted, no force potential is used.

For 3 or more templates, use template_extra:

    template_extra = '[1ATP][A][X][5],[8U18][B][Y][4]'

- Format: [PDB_CODE][PREDICTION_IDS][TEMPLATE_IDS][THRESHOLD]
- Extras support PDB codes only (no local file paths).
- Threshold is optional; if present, a force potential is applied.

### Affinity (Optional)

- `compute_affinity` — Whether to compute affinity predictions. Default: `true`. Set `false` to skip affinity calculations entirely.
- `affinity_ligand_id` — Ligand chain ID to use for affinity calculations. Defaults to `lig1_chain_id` if not set (only when it has a single chain ID; multi-copy ligands must set this explicitly, otherwise the pipeline errors).
- `affinity_mw_correction` — Whether to add the molecular weight correction to the affinity value head. Default: `false`.
- `diffusion_samples_affinity` — Number of diffusion samples for affinity prediction. Default: `5`.
- `sampling_steps_affinity` — Number of sampling steps for affinity prediction. Default: `200`.

### MSA Generation & Caching

Protein MSAs are generated with ColabFold (see [Setup](#3-colabfold-sequence-databases-for-msa-generation))
and cached by content hash so they're reused automatically across runs/ligand screens sharing
the same protein sequence(s):

- `msa_cache_dir`   — Cache directory (default: `./msa_cache`)
- `use_msa_cache`   — Set `false` to always regenerate MSAs instead of reusing the cache

### Boltz-2 Prediction Parameters

- `diffusion_samples` — Number of diffusion samples to generate. Default: `1` (Boltz-2 default).
- `seed` — Random seed(s) for reproducible predictions. A single integer (e.g. `42`), a comma-separated list (e.g. `'1,2,3'`, each generating a separate prediction with suffix `_s<N>`), or `null` for non-deterministic runs.
- `use_potentials` — Whether to use potentials to improve the physical quality of predictions (Boltz-2x). Default: `true`.
- `method` — Method conditioning to use. One of `afdb`, `boltz-1`, `md`, `electron microscopy`, `solution-nmr`, `x-ray diffraction`, `other`. Default: `null`.
- `recycling_steps` — Number of recycling steps. Default: `3` (Boltz-2 default).
- `sampling_steps` — Number of sampling steps. Default: `200` (Boltz-2 default).
- `step_scale` — Temperature-like scale for the diffusion sampling process; lower values increase diversity among samples (recommended range 1-2). Default: `1.638` (Boltz-2 default).
- `output_format` — Output structure format, `pdb` (default) or `mmcif`.
- `enable_visualisation` — Whether to generate ChimeraX visualizations. Default: `true`.
- `extra_config` — Additional Boltz-2 CLI arguments not otherwise exposed, e.g. `'--msa_pairing_strategy complete'`.

### Low-Memory Options

Only take effect when `use_lowmem_container = true` (uses the low-memory Boltz-2 fork built via `boltz2-lowmem.def`, see [Setup](#1-build-the-containers)):

- `use_lowmem_container` — Whether to use the low-memory Boltz-2 container. Default: `false`.
- `chunk_size_transition_z` — Transition Z chunk size. Default: `32`.
- `chunk_size_transition_msa` — Transition MSA chunk size. Default: `32`.
- `chunk_size_tri_attn` — Triangle attention chunk size. Default: `64` (upstream default: `128`).
- `triangle_mult_gate_nchunks` — Number of chunks for triangle multiplication. Default: `4`.
- `chunk_size_outer_product` — Outer product chunk size. Default: `4`.
- `chunk_size_threshold` — Maximum size before chunking is applied. Default: `384`.
- `use_bfloat16` — Whether to use bfloat16 precision. Default: `false`.

### Compute Resources

- `boltz2_gpu_model`, `boltz2_cpus`, `boltz2_memory` — Resources requested for the `RunBoltz` process. `boltz2_gpu_model` (`A30` or `A100`) and `gpu_queue` only take effect under a SLURM profile (e.g. `milton`); ignored for local execution.
- `colabfold_gpu_model`, `colabfold_cpus`, `colabfold_memory` — Resources requested for the `RunColabFold` process (MSA generation).
- Container paths (`boltz2_container`, `boltz2_lowmem_container`, `chimerax_container`, `colabfold_container`, `python_tools_container`) and model/database directories (`boltz2_models`, `colabfold_dbs`) — see [Setup](#setup).

## YAML Input Mode

You can supply a complete, prepared input YAML using --input_yaml.
> In this case, *all other sequence/ligand/constraint/template parameters are ignored*.

For YAML format details, see Boltz documentation:
https://github.com/jwohlwend/boltz/blob/main/docs/prediction.md

## Examples

> All examples below assume an execution profile is appended, e.g. `-profile standalone` or `-profile milton`.

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

## Citations

If you use this workflow, please cite it and the underlying tools it depends on:

1. Hardy, J.M. nf-boltz-2. *WorkflowHub* (2026). https://doi.org/10.48546/workflowhub.workflow.2318.1
2. Passaro, S., Corso, G., Wohlwend, J., Reveiz, M., Thaler, S., Somnath, V.R., Getz, N., Portnoi, T., Roy, J., Stark, H., Kwabi-Addo, D., Beaini, D., Jaakkola, T., Barzilay, R. Boltz-2: Towards Accurate and Efficient Binding Affinity Prediction. *bioRxiv* (2025). https://doi.org/10.1101/2025.06.14.659707
3. Litfin, T., Caley, J.S., Michie, K.A. LMI4Boltz: Optimising VRAM utilisation to predict large macromolecular complexes with consumer grade hardware. *bioRxiv* (2025). https://doi.org/10.1101/2025.10.29.684571
4. Di Tommaso, P., Chatzou, M., Floden, E.W., Barja, P.P., Palumbo, E., Notredame, C. Nextflow enables reproducible computational workflows. *Nature Biotechnology* 35, 316–319 (2017). https://doi.org/10.1038/nbt.3820
5. Kurtzer, G.M., Sochat, V., Bauer, M.W. Singularity: Scientific containers for mobility of compute. *PLOS ONE* 12, e0177459 (2017). https://doi.org/10.1371/journal.pone.0177459
6. Mirdita, M., Schütze, K., Moriwaki, Y., Heo, L., Ovchinnikov, S., Steinegger, M. ColabFold: making protein folding accessible to all. *Nature Methods* 19, 679–682 (2022). https://doi.org/10.1038/s41592-022-01488-1
7. Steinegger, M., Söding, J. MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. *Nature Biotechnology* 35, 1026–1028 (2017). https://doi.org/10.1038/nbt.3988
8. Pettersen, E.F., Goddard, T.D., Huang, C.C., Meng, E.C., Couch, G.S., Croll, T.I., Morris, J.H., Ferrin, T.E. UCSF ChimeraX: Structure visualization for researchers, educators, and developers. *Protein Science* 30, 70–82 (2021). https://doi.org/10.1002/pro.3943

## Credits

**nf-boltz-2** was originally written by [Joshua Hardy](https://github.com/JoshuaMHardy).

The A3M to CSV for Boltz script and ColabFold container definition were modified from work by Josh Caley (@jscgh) and the [nf-proteinfold](https://github.com/nf-core/proteinfold) project.
