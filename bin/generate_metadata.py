import json
import csv
import argparse
import numpy as np
import matplotlib.pyplot as plt
import os
import re
import pickle
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle
from matplotlib.collections import LineCollection
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import MaxNLocator
from glob import glob
import yaml as yaml_lib


def load_residue_offsets(directory):
    """Load residue_offsets from the YAML file in the results directory.
    Returns a dict mapping chain_id -> start_index, defaulting to empty dict.
    """
    yaml_files = glob(os.path.join(directory, '*.yaml'))
    for yaml_file in yaml_files:
        try:
            with open(yaml_file, 'r') as f:
                data = yaml_lib.safe_load(f)
            if data and 'residue_offsets' in data:
                return data['residue_offsets']
        except Exception:
            continue
    return {}


def detect_runname(directory):
    """Auto-detect RUNNAME from files in directory"""
    plddt_files = glob(os.path.join(directory, 'plddt_*_model_*.npz'))
    if not plddt_files:
        raise ValueError("No pLDDT files found in directory")
    base = os.path.basename(plddt_files[0])
    return base.split('_model_')[0].split('plddt_')[1]


def load_chains_data(directory, runname):
    """Load chain and residue information from npz file"""
    chains_file = os.path.join(directory, f"{runname}.npz")
    data = np.load(chains_file, allow_pickle=True)
    return data['chains'], data['residues']


def get_molecule_type_info():
    """Return molecule type mapping and labeling"""
    return {
        0: {'label': 'PROT', 'axis_label': 'Residue'},
        1: {'label': 'DNA', 'axis_label': 'Base'},
        2: {'label': 'RNA', 'axis_label': 'Base'},
        3: {'label': 'LIG', 'axis_label': 'Atom'}
    }


def create_residue_mapping(chains, residues):
    """Create mapping between atom-level indices and residue-level indices"""
    atom_to_residue_map = []
    residue_boundaries = []  # For tracking which atoms belong to which residue
    current_atom_idx = 0
    current_residue_idx = 0
    
    mol_type_info = get_molecule_type_info()
    
    for chain in chains:
        mol_type = chain['mol_type'].item()
        
        if mol_type == 3:  # Ligand - keep as atoms
            length = chain['atom_num'].item()
            for _ in range(length):
                atom_to_residue_map.append(current_residue_idx)
                residue_boundaries.append((current_atom_idx, current_atom_idx + 1))
                current_atom_idx += 1
                current_residue_idx += 1
        else:  # Protein/DNA/RNA - map atoms to residues
            start_res_idx = chain['res_idx'].item()
            end_res_idx = start_res_idx + chain['res_num'].item()
            
            for res_idx in range(start_res_idx, end_res_idx):
                if res_idx < len(residues):
                    residue = residues[res_idx]
                    if residue['is_standard']:
                        # Standard residue: one atom
                        atom_to_residue_map.append(current_residue_idx)
                        residue_boundaries.append((current_atom_idx, current_atom_idx + 1))
                        current_atom_idx += 1
                    else:
                        # Non-standard residue: multiple atoms map to one residue
                        atom_count = residue['atom_num']
                        for _ in range(atom_count):
                            atom_to_residue_map.append(current_residue_idx)
                        residue_boundaries.append((current_atom_idx, current_atom_idx + atom_count))
                        current_atom_idx += atom_count
                    current_residue_idx += 1
    
    return atom_to_residue_map, residue_boundaries, current_residue_idx


def average_data_by_residue(data, residue_boundaries):
    """Average any data type (pLDDT, PAE, PDE) by residue using residue boundaries"""
    if data.ndim == 1:
        # 1D data (like pLDDT)
        averaged_data = []
        for start, end in residue_boundaries:
            if end > start:
                averaged_data.append(np.mean(data[start:end]))
            else:
                averaged_data.append(0.0)  # Fallback for empty ranges
        return np.array(averaged_data)
    
    elif data.ndim == 2:
        # 2D data (like PAE, PDE matrices)
        n_residues = len(residue_boundaries)
        averaged_matrix = np.zeros((n_residues, n_residues))
        
        for i, (start_i, end_i) in enumerate(residue_boundaries):
            for j, (start_j, end_j) in enumerate(residue_boundaries):
                if end_i > start_i and end_j > start_j:
                    # Average the submatrix
                    submatrix = data[start_i:end_i, start_j:end_j]
                    averaged_matrix[i, j] = np.mean(submatrix)
        
        return averaged_matrix
    
    else:
        raise ValueError(f"Unsupported data dimensionality: {data.ndim}")


def get_processed_chain_ranges(chains, residues, residue_offsets=None):
    """Get chain ranges after residue-level processing.
    If residue_offsets is provided, includes start_index for each chain.
    Returns: (chain_ranges, total_idx) where chain_ranges is list of
    (chain_id, mol_label, start_idx, end_idx, start_index) tuples.
    start_index is the user's residue numbering offset (default 1).
    """
    if residue_offsets is None:
        residue_offsets = {}
    chain_ranges = []
    current_idx = 0
    mol_type_info = get_molecule_type_info()
    
    for chain in chains:
        chain_id = chain['name'].item()
        mol_type = chain['mol_type'].item()
        mol_label = mol_type_info[mol_type]['label']
        chain_start_idx = current_idx
        start_index = residue_offsets.get(chain_id, 1)
        
        if mol_type == 3:  # Ligand - count atoms
            length = chain['atom_num'].item()
            current_idx += length
        else:  # Protein/DNA/RNA - count residues
            start_res_idx = chain['res_idx'].item()
            end_res_idx = start_res_idx + chain['res_num'].item()
            current_idx += (end_res_idx - start_res_idx)
        
        chain_ranges.append((chain_id, mol_label, chain_start_idx, current_idx, start_index))
    
    return chain_ranges, current_idx


def create_confidence_legend():
    """Create confidence legend elements"""
    return [
        Rectangle((0, 0), 1, 1, facecolor=(0.051, 0.341, 0.827), alpha=0.5, label='Very high (pLDDT > 90)'),
        Rectangle((0, 0), 1, 1, facecolor=(0.416, 0.796, 0.945), alpha=0.5, label='Confident (90 > pLDDT > 70)'),
        Rectangle((0, 0), 1, 1, facecolor=(0.996, 0.851, 0.212), alpha=0.5, label='Low (70 > pLDDT > 50)'),
        Rectangle((0, 0), 1, 1, facecolor=(0.992, 0.490, 0.302), alpha=0.5, label='Very low (pLDDT < 50)')
    ]


def add_confidence_background(ax):
    """Add confidence background bands"""
    ax.axhspan(0, 50, color=(0.992, 0.490, 0.302), alpha=0.3, zorder=0)
    ax.axhspan(50, 70, color=(0.996, 0.851, 0.212), alpha=0.3, zorder=0)
    ax.axhspan(70, 90, color=(0.416, 0.796, 0.945), alpha=0.3, zorder=0)
    ax.axhspan(90, 100, color=(0.051, 0.341, 0.827), alpha=0.3, zorder=0)


def plot_plddt_chain(ax, pl_values, chain_id, mol_label, is_first_chain=False, start_index=1):
    """Plot pLDDT data for a single chain"""

    mol_type_info = get_molecule_type_info()
    x_vals = np.arange(start_index, start_index + len(pl_values))
    
    add_confidence_background(ax)
    ax.xaxis.set_major_locator(MaxNLocator(steps = [1, 2, 5, 10]))

    # Handle single data point differently
    if len(pl_values) == 1:
        ax.scatter(x_vals[0], pl_values[0], s=10, zorder=5, color=(0.2, 0.2, 0.2))
        ax.set_xlim(x_vals[0] - 0.5, x_vals[0] + 0.5)
        ax.set_xticks([x_vals[0]])
    else:
        # Multiple points: use LineCollection method
        points = np.array([x_vals, pl_values]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        lc = LineCollection(
            segments,
            color=(0.2,0.2,0.2),
            linewidths=2)
        ax.add_collection(lc)
        
        ax.set_xlim(x_vals[0], x_vals[-1])
  
    ax.set_ylim(0, 100)
    ax.yaxis.set_major_locator(MultipleLocator(10))
    
    if is_first_chain:
        ax.set_ylabel("pLDDT Confidence")
    else:
        ax.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)
    
    # Get axis label from molecule type
    xlabel = next(info['axis_label'] for key, info in mol_type_info.items() 
                 if info['label'] == mol_label)
    ax.set_xlabel(xlabel)
    ax.set_title(f"Chain {chain_id}\n({mol_label})", fontsize=12)


def plot_single_model_plddt(directory, runname, output_dir, chains, residues, model_num, residue_offsets=None):
    """Multi-panel pLDDT plot with horizontal layout and averaged non-standard residues"""
    file = os.path.join(directory, f'plddt_{runname}_model_{model_num}.npz')
    data = np.load(file)
    raw_plddt = data['plddt'] * 100
    
    # Process pLDDT data
    atom_to_residue_map, residue_boundaries, total_residues = create_residue_mapping(chains, residues)
    processed_plddt = average_data_by_residue(raw_plddt, residue_boundaries)
    chain_ranges, _ = get_processed_chain_ranges(chains, residues, residue_offsets)
    
    n_chains = len(chain_ranges)
    chain_lengths = [max(end - start, 1) for _, _, start, end, _ in chain_ranges]
    width_ratios = [max(length, 0.3) for length in chain_lengths]
    
    fig, axes = plt.subplots(1, n_chains, figsize=(16, 5), sharey=True, 
                            gridspec_kw={'width_ratios': width_ratios})
    if n_chains == 1:
        axes = [axes]
    
    for i, (chain_id, mol_label, start, end, start_index) in enumerate(chain_ranges):
        pl_values = processed_plddt[start:end]
        plot_plddt_chain(axes[i], pl_values, chain_id, mol_label, is_first_chain=(i == 0), start_index=start_index)
    
    # Add legend
    legend_elements = create_confidence_legend()
    fig.legend(
        handles=legend_elements,
        title='pLDDT confidence (non-standard residues averaged):',
        loc='lower center', 
        bbox_to_anchor=(0.5, 0.01), 
        ncol=4, 
        fontsize=10, 
        framealpha=0.9, 
        title_fontsize=11
    )
    
    plt.suptitle(f"pLDDT Confidence - Model {model_num}", 
                fontsize=15, fontweight='bold')
    
    plt.subplots_adjust(
        left=0.06, right=0.95, top=0.85, bottom=0.25, wspace=0.05
    )
    
    plt.savefig(os.path.join(output_dir, f'{runname}_model_{model_num}_plddt.png'),
                dpi=150, bbox_inches='tight', pad_inches=0.1)
    plt.close()


def plot_plddt_chain_multi_model(ax, all_model_values, chain_id, mol_label, is_first_chain=False, start_index=1):
    """Plot pLDDT data for multiple models on a single chain"""
    mol_type_info = get_molecule_type_info()
    
    add_confidence_background(ax)

    ax.xaxis.set_major_locator(MaxNLocator(steps = [1, 2, 5, 10]))

    # Handle single data point differently
    if len(all_model_values[0]) == 1:
        # For single points, use scatter plots with different colors per model
        colors = plt.cm.tab10(np.linspace(0, 1, len(all_model_values)))
        for model_idx, pl_values in enumerate(all_model_values):
            ax.scatter(start_index, pl_values[0], color=colors[model_idx], s=10, 
                      label=f'Model {model_idx}', zorder=5, alpha=0.8)
        ax.set_xlim(start_index - 0.5, start_index + 0.5)
        ax.set_xticks([start_index])
    else:
        # Multiple points: use regular line plots with different colors per model
        x_vals = np.arange(start_index, start_index + len(all_model_values[0]))
        colors = plt.cm.tab10(np.linspace(0, 1, len(all_model_values)))
        
        for model_idx, pl_values in enumerate(all_model_values):
            ax.plot(x_vals, pl_values, color=colors[model_idx], alpha=0.8, 
                   label=f'Model {model_idx}', linewidth=2)
        
        ax.set_xlim(x_vals[0], x_vals[-1])

    ax.set_ylim(0, 100)
    ax.yaxis.set_major_locator(MultipleLocator(10))
    
    if is_first_chain:
        ax.set_ylabel("pLDDT Confidence")
        # Add model legend only to the first subplot
        ax.legend(loc='upper left', fontsize=9)
    else:
        ax.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)
    
    # Get axis label from molecule type
    xlabel = next(info['axis_label'] for key, info in mol_type_info.items() 
                 if info['label'] == mol_label)
    ax.set_xlabel(xlabel)
    ax.set_title(f"Chain {chain_id}\n({mol_label})", fontsize=12)


def plot_all_models_plddt(directory, runname, output_dir, chains, residues, residue_offsets=None):
    """Plot pLDDT values for all models with horizontal layout and averaged non-standard residues"""
    plddt_files = sorted(glob(os.path.join(directory, f'plddt_{runname}_model_*.npz')))
    
    # Process data
    atom_to_residue_map, residue_boundaries, total_residues = create_residue_mapping(chains, residues)
    chain_ranges, _ = get_processed_chain_ranges(chains, residues, residue_offsets)
    
    n_chains = len(chain_ranges)
    n_models = len(plddt_files)
    
    # Load and process all model data
    all_processed_plddt = []
    for file in plddt_files:
        data = np.load(file)
        raw_plddt = data['plddt'] * 100
        processed_plddt = average_data_by_residue(raw_plddt, residue_boundaries)
        all_processed_plddt.append(processed_plddt)
    
    # Calculate width ratios based on chain lengths (same as single model plot)
    chain_lengths = [max(end - start, 1) for _, _, start, end, _ in chain_ranges]
    width_ratios = [max(length, 0.3) for length in chain_lengths]
    
    # Create horizontal subplots with proportional widths (matching single model layout)
    fig, axes = plt.subplots(1, n_chains, figsize=(16, 5), sharey=True, 
                            gridspec_kw={'width_ratios': width_ratios})
    if n_chains == 1:
        axes = [axes]
    
    # Plot each chain
    for i, (chain_id, mol_label, start, end, start_index) in enumerate(chain_ranges):
        # Get pLDDT values for all models for this chain
        all_model_chain_values = []
        for processed_plddt in all_processed_plddt:
            pl_values = processed_plddt[start:end]
            all_model_chain_values.append(pl_values)
        
        plot_plddt_chain_multi_model(axes[i], all_model_chain_values, chain_id, mol_label, 
                                    is_first_chain=(i == 0), start_index=start_index)
    
    # Add confidence legend (same as single model plot)
    legend_elements = create_confidence_legend()
    fig.legend(
        handles=legend_elements,
        title='pLDDT confidence (non-standard residues averaged):',
        loc='lower center', 
        bbox_to_anchor=(0.5, 0.01), 
        ncol=4, 
        fontsize=10, 
        framealpha=0.9, 
        title_fontsize=11
    )
    
    plt.suptitle("pLDDT Confidence Scores - All Models", 
                fontsize=15, fontweight='bold')
    
    # Use same layout adjustments as single model plot
    plt.subplots_adjust(
        left=0.06, right=0.95, top=0.85, bottom=0.25, wspace=0.05
    )
    
    plt.savefig(os.path.join(output_dir, f'{runname}_all_models_plddt.png'), 
               dpi=150, bbox_inches='tight', pad_inches=0.1)
    plt.close()


def plot_2d_matrix(directory, runname, output_dir, matrix_type, model_num, chains, residues, residue_offsets=None):
    """Plot PAE or PDE matrix with non-standard residue averaging and optimized layout"""
    file = os.path.join(directory, f"{matrix_type}_{runname}_model_{model_num}.npz")
    if not os.path.exists(file):
        return
        
    data = np.load(file)
    raw_matrix = data[matrix_type]
    
    # Process matrix data
    atom_to_residue_map, residue_boundaries, total_residues = create_residue_mapping(chains, residues)
    processed_matrix = average_data_by_residue(raw_matrix, residue_boundaries)
    chain_ranges, _ = get_processed_chain_ranges(chains, residues, residue_offsets)
    
    # Create figure with optimized size - made it square to reduce whitespace
    fig, ax = plt.subplots(figsize=(8, 8))
    
    if matrix_type == 'pae':
        im = ax.imshow(processed_matrix, cmap='GnBu_r', origin='lower', vmin=0, vmax=30, 
                      extent=[-0.5, total_residues-0.5, -0.5, total_residues-0.5])
        cbar_label = 'Expected position error (Å)'
        explanation = (
            f"{matrix_type.upper()} measures confidence in relative positions. "
            "Dark blue = high confidence, light colors = uncertainty. Non-standard residues averaged."
        )
    else:  # PDE
        im = ax.imshow(processed_matrix, cmap='GnBu_r', origin='lower',
                      extent=[-0.5, total_residues-0.5, -0.5, total_residues-0.5])
        cbar_label = 'Expected distance error (Å)'
        explanation = (
            f"{matrix_type.upper()} measures distance confidence. "
            "Dark blue = high confidence, light colors = uncertainty. Non-standard residues averaged."
        )
    
    # Add chain boundary lines
    for i, (chain_id, mol_label, start, end, start_index) in enumerate(chain_ranges):
        if i > 0:
            ax.axhline(start - 0.5, color='red', linestyle='-', alpha=0.8, linewidth=1.5)
            ax.axvline(start - 0.5, color='red', linestyle='-', alpha=0.8, linewidth=1.5)
    
    # Set plot limits
    ax.set_xlim(-0.5, total_residues - 0.5)
    ax.set_ylim(-0.5, total_residues - 0.5)
    
    # Add axis labels
    ax.set_xlabel("Scored residue", fontsize=12)
    ax.set_ylabel("Aligned residue", fontsize=12)
    
    # Reduced label offset to minimize whitespace
    label_offset = total_residues * 0.02
    
    for chain_id, mol_label, start, end, start_index in chain_ranges:
        mid_pos = (start + end) / 2
        # Show start_index in label if non-default
        offset_info = f"\n(from {start_index})" if start_index != 1 else ""
        
        # Top labels - closer to plot
        ax.text(mid_pos, total_residues + label_offset, f"{chain_id}\n{mol_label}{offset_info}", 
                ha='center', va='bottom', fontsize=9,
                bbox=dict(boxstyle="round,pad=0.2", facecolor="white", 
                         alpha=0.9, edgecolor='gray'))
        
        # Right labels - closer to plot
        ax.text(total_residues + label_offset, mid_pos, f"{chain_id} - {mol_label}{offset_info}",
                ha='left', va='center', fontsize=9,
                bbox=dict(boxstyle="round,pad=0.2", facecolor="white", 
                         alpha=0.9, edgecolor='gray'))
    
    # Add colorbar with better positioning
    cbar = fig.colorbar(im, ax=ax, orientation='horizontal', 
                       shrink=0.8, pad=0.1, aspect=35, fraction=0.046)
    cbar.set_label(cbar_label, fontsize=11)
    cbar.ax.tick_params(labelsize=10)
    
    # Position title at the very top, above chain labels
    fig.suptitle(f"{matrix_type.upper()} - Model {model_num}", 
                fontsize=14, fontweight='bold', y=0.98)
    
    # Position explanation text more compactly
    fig.text(0.5, 0.03, explanation, ha='center', va='bottom', fontsize=9,
             bbox=dict(boxstyle="round,pad=0.3", facecolor="lightgray", alpha=0.8),
             transform=fig.transFigure, wrap=True)
    
    # Adjusted layout to accommodate title at the top
    plt.subplots_adjust(
        top=0.9,       # Reduced top margin to make room for title above chain labels
        bottom=0.12     # Reduced bottom margin
    )
    
    plt.savefig(os.path.join(output_dir, f'{runname}_model_{model_num}_{matrix_type}.png'), 
               dpi=150, bbox_inches='tight', pad_inches=0.1)
    plt.close()


def plot_iptm_matrix(directory, runname, output_dir, model_num, chain_ranges, iptm_matrix):
    """Plot chain interaction ipTM matrix"""
    plt.figure(figsize=(8, 8))
    chain_ids = [f"{chain_id}\n{mol_label}" for chain_id, mol_label, _, _, _ in chain_ranges]
    num_chains = len(chain_ids)
    symmetric_matrix = np.maximum(iptm_matrix, iptm_matrix.T)
    masked_matrix = np.ma.masked_where(symmetric_matrix == 0, symmetric_matrix)
    cmap = plt.get_cmap('viridis')
    im = plt.imshow(masked_matrix, cmap=cmap, origin='lower', vmin=0, vmax=1, aspect='equal')
    plt.xticks(range(num_chains), chain_ids, fontsize=10)
    plt.yticks(range(num_chains), chain_ids, fontsize=10)
    cbar = plt.colorbar(im, label='Interface Confidence (ipTM)', shrink=0.8, pad=0.12, fraction=0.05, aspect=35, orientation='horizontal')
    cbar.set_ticks(np.arange(0, 1.1, 0.2))
    
    for i in range(num_chains):
        for j in range(num_chains):
            if symmetric_matrix[i, j] > 0:
                plt.text(j, i, f'{symmetric_matrix[i, j]:.2f}', 
                        ha='center', va='center', 
                        color='white' if symmetric_matrix[i, j] < 0.5 else 'black')
            else:
                plt.text(j, i, 'N/A', ha='center', va='center', color='grey', style='italic')
    
    plt.title(f"Interface Confidence - Model {model_num}")
    plt.xlabel("Interacting Chain", fontsize=12)
    plt.ylabel("Interacting Chain", fontsize=12)
    
    explanation = (
        "ipTM scores estimate interface prediction confidence (0-1).\n"
        "Values show mutual interface quality between chain pairs.\n"
        "Higher scores indicate more reliable interactions."
    )
    plt.figtext(0.5, 0.01, explanation, ha='center', va='bottom', fontsize=9,
                bbox=dict(boxstyle="round,pad=0.5", facecolor="lightgray", alpha=0.8))
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.12)
    plt.savefig(os.path.join(output_dir, f'{runname}_model_{model_num}_iptm.png'), dpi=150, bbox_inches='tight')
    plt.close()


def process_confidence_data(directory, runname, model_num, chains, residue_offsets=None):
    """Process confidence JSON file and return metadata and matrix with ordered keys"""
    json_file = os.path.join(directory, f'confidence_{runname}_model_{model_num}.json')
    if not os.path.exists(json_file):
        return None, None
    
    with open(json_file) as f:
        data = json.load(f)
    
    # Create chain ID mapping
    chain_id_map = {}
    chain_ranges, _ = get_processed_chain_ranges(chains, [], residue_offsets)
    
    for idx, (chain_id, mol_label, start, end, _start_index) in enumerate(chain_ranges):
        chain_data = next(c for c in chains if c['name'].item() == chain_id)
        chain_id_map[chain_data['asym_id'].item()] = {
            'id': chain_id,
            'type': mol_label,
            'index': idx
        }
    
    # Define ordered fieldnames following JSON structure
    ordered_fieldnames = [
        'pred_name',
        'model',
        'confidence_score',
        'ptm',
        'iptm',
        'ligand_iptm',
        'protein_iptm',
        'complex_plddt',
        'complex_iplddt',
        'complex_pde',
        'complex_ipde'
    ]
    
    # Add chains_ptm fields in alphabetical chain order
    chain_ids = sorted([info['id'] for info in chain_id_map.values()])
    for chain_id in chain_ids:
        ordered_fieldnames.append(f'ptm_{chain_id}')
    
    # Add pair_chains_iptm fields in alphabetical pair order
    for src_chain in chain_ids:
        for dst_chain in chain_ids:
            ordered_fieldnames.append(f'iptm_{src_chain}_{dst_chain}')
    
    # Process metadata in order
    metadata = {}
    
    # Main metrics (in order)
    metadata['pred_name'] = runname
    metadata['model'] = model_num
    metadata['confidence_score'] = round(data.get('confidence_score', 0), 3)
    metadata['ptm'] = round(data.get('ptm', 0), 3)
    metadata['iptm'] = round(data.get('iptm', 0), 3)
    metadata['ligand_iptm'] = round(data.get('ligand_iptm', 0), 3)
    metadata['protein_iptm'] = round(data.get('protein_iptm', 0), 3)
    metadata['complex_plddt'] = round(data.get('complex_plddt', 0), 3)
    metadata['complex_iplddt'] = round(data.get('complex_iplddt', 0), 3)
    metadata['complex_pde'] = round(data.get('complex_pde', 0), 3)
    metadata['complex_ipde'] = round(data.get('complex_ipde', 0), 3)
    
    # Process chains_ptm in chain order
    chains_ptm = data.get('chains_ptm', {})
    for eid, value in chains_ptm.items():
        chain_info = chain_id_map.get(int(eid))
        if chain_info and float(value) != 0.0:
            metadata[f'ptm_{chain_info["id"]}'] = round(float(value), 3)
    
    # Create ipTM matrix
    num_chains = len(chain_id_map)
    iptm_matrix = np.zeros((num_chains, num_chains), dtype=np.float32)
    
    # Process pair_chains_iptm in chain order
    pair_iptm = data.get('pair_chains_iptm', {})
    for src_eid, targets in pair_iptm.items():
        src_chain = chain_id_map.get(int(src_eid))
        if not src_chain:
            continue
            
        for dst_eid, value in targets.items():
            dst_chain = chain_id_map.get(int(dst_eid))
            value = float(value)
            if value == 0.0 or not dst_chain:
                continue
                
            metadata[f'iptm_{src_chain["id"]}_{dst_chain["id"]}'] = round(value, 3)
            iptm_matrix[src_chain['index'], dst_chain['index']] = value
    
    # Store the ordered fieldnames for later use
    metadata['_ordered_fieldnames'] = ordered_fieldnames
    
    return metadata, iptm_matrix

# Utility functions for affinity data (unchanged)
def affinity_ensemble_row(model_idx, pred_log_ic50, prob_binder):
    """Convert raw log(IC50) and probability to multiple human-readable columns"""
    ic50_um = 10 ** pred_log_ic50
    pIC50 = 6 - pred_log_ic50
    deltaG = (6 - pred_log_ic50) * 1.364
    return {
        "ensemble": model_idx,
        "prob_binder": f"{prob_binder:.3f}",
        "pred_IC50_μM": f"{ic50_um:.3g}",
        "pred_pIC50": f"{pIC50:.2f}",
        "pred_deltaG": f"{deltaG:.2f}",
    }


def read_affinity_ensembles(directory, runname):
    """Extract all affinity ensemble predictions"""
    affinity_json = os.path.join(directory, f'affinity_{runname}.json')
    if not os.path.exists(affinity_json):
        return []
    
    with open(affinity_json) as f:
        affinity_data = json.load(f)
    
    ensemble_rows = []
    idx = 0
    while True:
        suffix = "" if idx == 0 else str(idx)
        pred_key = f"affinity_pred_value{suffix}"
        prob_key = f"affinity_probability_binary{suffix}"
        if pred_key in affinity_data and prob_key in affinity_data:
            pred = float(affinity_data[pred_key])
            prob = float(affinity_data[prob_key])
            ensemble_rows.append(affinity_ensemble_row(idx, pred, prob))
            idx += 1
        else:
            break
    return ensemble_rows

def write_pickle_files(directory, runname, output_dir):
    """Write pLDDT, PAE, and PDE data to pickle files"""

    # Process each data type
    for matrix_type in ['plddt', 'pae', 'pde']:
        # Find all npz files of this type
        npz_files = sorted([f for f in os.listdir(directory) 
                           if f.startswith(f"{matrix_type}_{runname}") and f.endswith(".npz")])
        
        for file_name in npz_files:
            # Extract model number from filename
            match = re.search(rf"{matrix_type}_{runname}_model_(\d+)\.npz", file_name)
            if not match:
                continue
            model_num = int(match.group(1))
            
            file_path = os.path.join(directory, file_name)
            if not os.path.exists(file_path):
                continue
                
            # Load and process data
            data = np.load(file_path)
            data = data[matrix_type]

            # Write to pickle file in ChimeraX-compatible format
            pkl_filename = f"{runname}_{matrix_type}_model_{model_num}.pkl"
            pkl_filepath = os.path.join(output_dir, pkl_filename)
            
            # **KEY FIX**: Wrap in dictionary with expected key for ChimeraX compatibility
            if matrix_type in ['pae', 'pde']:
                # For PAE and PDE matrices, use 'predicted_aligned_error' key
                pickle_data = {'predicted_aligned_error': data}
            else:
                # For pLDDT (1D data), you might want a different key or structure
                # ChimeraX may not have specific support for pLDDT pickle files
                pickle_data = {'plddt': data}
            
            with open(pkl_filepath, 'wb') as pkl_file:
                pickle.dump(pickle_data, pkl_file)
            
            print(f"Saved {matrix_type} data to {pkl_filepath}")
    
    print(f"All pickle files saved to {output_dir}")

def main():
    parser = argparse.ArgumentParser(description='Generate Boltz-2 analysis plots')
    parser.add_argument('input_dir', help='Path to directory containing Boltz-2 results')
    parser.add_argument('-o', '--output', default='plots', help='Output directory')
    parser.add_argument('--write-pkl', action='store_true', 
                       help='Write data to pickle files for external analysis')
    args = parser.parse_args()
    
    os.makedirs(args.output, exist_ok=True)
    runname = detect_runname(args.input_dir)
    chains, residues = load_chains_data(args.input_dir, runname)
    
    # Load residue offsets from YAML if available
    residue_offsets = load_residue_offsets(args.input_dir)
    if residue_offsets:
        print(f"Loaded residue offsets: {residue_offsets}")
    
    # Write  pickle files if requested
    if args.write_pkl:
        print("Writing pickle files...")
        write_pickle_files(args.input_dir, runname, args.output)
    
    # Get processed chain ranges for metadata
    chain_ranges, _ = get_processed_chain_ranges(chains, residues, residue_offsets)
    
    # CSV metadata collection
    all_metadata = []
    ordered_fieldnames = None
    
    # Plot summary pLDDT plot
    plot_all_models_plddt(args.input_dir, runname, args.output, chains, residues, residue_offsets)
    
    # Find number of models
    num_models = len(glob(os.path.join(args.input_dir, f'plddt_{runname}_model_*.npz')))
    
    for model_num in range(num_models):
        # Per-model pLDDT plot
        plot_single_model_plddt(args.input_dir, runname, args.output, chains, residues, model_num, residue_offsets)
        
        # PAE and PDE plots with averaging
        plot_2d_matrix(args.input_dir, runname, args.output, 'pae', model_num, chains, residues, residue_offsets)
        plot_2d_matrix(args.input_dir, runname, args.output, 'pde', model_num, chains, residues, residue_offsets)
        
        # Process confidence data and ipTM plot
        metadata, iptm_matrix = process_confidence_data(args.input_dir, runname, model_num, chains, residue_offsets)
        if metadata and iptm_matrix is not None:
            # Extract ordered fieldnames from first metadata entry
            if ordered_fieldnames is None:
                ordered_fieldnames = metadata.pop('_ordered_fieldnames', None)
            else:
                metadata.pop('_ordered_fieldnames', None)  # Remove from subsequent entries
            
            all_metadata.append(metadata)
            plot_iptm_matrix(args.input_dir, runname, args.output, model_num, chain_ranges, iptm_matrix)
    
    # Write CSV metadata with proper ordering
    if all_metadata:
        csv_path = os.path.join(args.output, f'{runname}_metadata.csv')
        
        if ordered_fieldnames:
            # Use the ordered fieldnames we defined
            fieldnames = ordered_fieldnames
        else:
            # Fallback to the old method if something goes wrong
            fieldnames = sorted(set(key for md in all_metadata for key in d.keys()))
        
        with open(csv_path, 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(all_metadata)
    
    # Write affinity CSV (unchanged)
    affinity_ensembles = read_affinity_ensembles(args.input_dir, runname)
    if affinity_ensembles:
        affinity_csv_path = os.path.join(args.output, f'{runname}_affinity.csv')
        fieldnames = ["ensemble", "pred_IC50_μM", "pred_pIC50", "pred_deltaG", "prob_binder"]
        with open(affinity_csv_path, 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(affinity_ensembles)

if __name__ == '__main__':
    main()
