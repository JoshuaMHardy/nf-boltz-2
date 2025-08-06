import json
import csv
import argparse
import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle
from matplotlib.collections import LineCollection
from glob import glob

def detect_runname(directory):
    """Auto-detect RUNNAME from files in directory"""
    plddt_files = glob(os.path.join(directory, 'plddt_*_model_*.npz'))
    if not plddt_files:
        raise ValueError("No pLDDT files found in directory")
    base = os.path.basename(plddt_files[0])
    return base.split('_model_')[0].split('plddt_')[1]

def affinity_ensemble_row(model_idx, pred_log_ic50, prob_binder):
    """Convert raw log(IC50) and probability to multiple human-readable columns for a single ensemble component."""
    # IC50 in μM (3 decimals), pIC50 (2 decimals), deltaG in kcal/mol (2 decimal), probability (3 decimals)
    ic50_um = 10 ** pred_log_ic50
    pIC50 = 6 - pred_log_ic50
    deltaG = (6 - pred_log_ic50) * 1.364
    return {
        "ensemble": model_idx,
        "pred_IC50_μM": f"{ic50_um:.3g}",       # e.g., "0.001", "0.1", "5.6", "120"
        "pred_pIC50": f"{pIC50:.2f}",           # e.g., "9.00"
        "pred_deltaG": f"{deltaG:.2f}",         # e.g., "8.18"
        "prob_binder": f"{prob_binder:.3f}",    # e.g., "0.126"
    }
def read_affinity_ensembles(directory, runname):
    """Extract all affinity ensemble predictions and convert for reporting."""
    affinity_json = os.path.join(directory, f'affinity_{runname}.json')
    if not os.path.exists(affinity_json):
        return []
    with open(affinity_json) as f:
        affinity_data = json.load(f)
    # Find all unique ensemble IDs by suffix
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

def load_chains_data(directory, runname):
    """Load chain and residue information from npz file"""
    chains_file = os.path.join(directory, f"{runname}.npz")
    data = np.load(chains_file, allow_pickle=True)
    return data['chains'], data['residues']

def parse_chain_ranges_with_residues(chains, residues):
    """Calculate residue ranges for each chain accounting for non-standard residues"""
    chain_ranges = []
    non_standard_ranges = []
    current_idx = 0
    mol_type_names = {0: 'PROT', 1: 'DNA', 2: 'RNA', 3: 'LIG'}
    
    for chain in chains:
        chain_id = chain['name'].item()
        mol_type = chain['mol_type'].item()
        mol_label = mol_type_names[mol_type]
        
        if mol_type == 3:  # Ligand
            length = chain['atom_num'].item()
            chain_ranges.append((chain_id, mol_label, current_idx, current_idx + length))
            current_idx += length
        else:  # Protein/DNA/RNA
            start_res_idx = chain['res_idx'].item()
            end_res_idx = start_res_idx + chain['res_num'].item()
            chain_start_idx = current_idx
            
            for res_idx in range(start_res_idx, end_res_idx):
                if res_idx < len(residues):
                    residue = residues[res_idx]
                    if residue['is_standard']:
                        current_idx += 1
                    else:
                        atom_count = residue['atom_num']
                        res_name = residue['name'].item()
                        non_standard_ranges.append((res_name, res_idx, current_idx, current_idx + atom_count))
                        current_idx += atom_count
            
            chain_ranges.append((chain_id, mol_label, chain_start_idx, current_idx))
    
    return chain_ranges, non_standard_ranges, current_idx

def create_alphafold_colormap():
    """Create corrected AlphaFold colormap"""
    colors = [
        (0.0, (0.992, 0.490, 0.302)),   # Red/orange
        (50/100, (0.992, 0.490, 0.302)),
        (50/100, (0.996, 0.851, 0.212)), # Yellow
        (70/100, (0.996, 0.851, 0.212)),
        (70/100, (0.416, 0.796, 0.945)), # Light blue
        (90/100, (0.416, 0.796, 0.945)),
        (90/100, (0.051, 0.341, 0.827)), # Dark blue
        (1.0, (0.051, 0.341, 0.827))
    ]
    
    return LinearSegmentedColormap.from_list(
        'alphafold_corrected',
        colors,
        N=256
    )

def add_confidence_legend():
    """Add confidence color legend to the plot"""
    legend_elements = [
        Rectangle((0, 0), 1, 1, facecolor=(0.051, 0.341, 0.827), label='Very high (pLDDT > 90)'),
        Rectangle((0, 0), 1, 1, facecolor=(0.416, 0.796, 0.945), label='Confident (90 > pLDDT > 70)'),
        Rectangle((0, 0), 1, 1, facecolor=(0.996, 0.851, 0.212), label='Low (70 > pLDDT > 50)'),
        Rectangle((0, 0), 1, 1, facecolor=(0.992, 0.490, 0.302), label='Very low (pLDDT < 50)')
    ]
    plt.legend(handles=legend_elements, title='Model Confidence:', 
              loc='upper right', bbox_to_anchor=(1, 1))

def add_confidence_background():
    """Add horizontal confidence bands to background"""
    plt.axhspan(0, 50, color=(0.992, 0.490, 0.302), alpha=0.1, zorder=0)  # Very low
    plt.axhspan(50, 70, color=(0.996, 0.851, 0.212), alpha=0.1, zorder=0)  # Low
    plt.axhspan(70, 90, color=(0.416, 0.796, 0.945), alpha=0.1, zorder=0)  # Confident
    plt.axhspan(90, 100, color=(0.051, 0.341, 0.827), alpha=0.1, zorder=0)  # Very high

def add_chain_labels_1d(chain_ranges, non_standard_ranges):
    """Add chain labels and non-standard residue markers to 1D plots"""
    # Chain boundaries - now dark grey
    for i, (chain_id, mol_label, start, end) in enumerate(chain_ranges):
        if i > 0:
            plt.axvline(start, color='#333333', linestyle='-', alpha=0.8, linewidth=2)
        
        # Chain labels
        plt.text((start + end)/2, 5, f"{chain_id}\n{mol_label}", 
                ha='center', va='bottom', fontsize=8, 
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
    
    # Non-standard residues - also dark grey
    for res_name, res_idx, start, end in non_standard_ranges:
        plt.axvline(start, color='#333333', linestyle='--', alpha=0.8, linewidth=1)
        plt.axvline(end, color='#333333', linestyle='--', alpha=0.8, linewidth=1)
        plt.text((start + end)/2, 15, f"{res_name}{res_idx}", 
                ha='center', va='bottom', fontsize=7,
                bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.8))

def plot_2d_matrix(directory, runname, output_dir, matrix_type, model_num, chain_ranges, non_standard_ranges, total_length):
    """Plot PAE or PDE matrix with compact layout and proper spacing"""
    file = os.path.join(directory, f"{matrix_type}_{runname}_model_{model_num}.npz")
    data = np.load(file)
    matrix = data[matrix_type][:total_length, :total_length]
    
    fig, ax = plt.subplots(figsize=(8, 8))  # More compact
    
    if matrix_type == 'pae':
        im = ax.imshow(matrix, cmap='GnBu_r', origin='lower', vmin=0, vmax=30, 
                      extent=[-0.5, total_length-0.5, -0.5, total_length-0.5])
        cbar_label = 'Expected prediction error (Å)'
        explanation = (
            "PAE shows expected positional error between residue pairs.\n"
            "Dark blue = low error (high confidence), light colors = higher uncertainty."
        )
    else:  # PDE
        matrix = np.log1p(matrix)
        im = ax.imshow(matrix, cmap='GnBu_r', origin='lower', 
                      extent=[-0.5, total_length-0.5, -0.5, total_length-0.5])
        cbar_label = 'log(Expected Distance Error (Å))'
        explanation = (
            "PDE estimates distance constraint reliability (log scale).\n"
            "Dark blue = confident predictions, light colors = uncertainty."
        )
    
    # Add chain boundary lines in RED
    for i, (chain_id, mol_label, start, end) in enumerate(chain_ranges):
        if i > 0:
            ax.axhline(start - 0.5, color='red', linestyle='-', alpha=0.8, linewidth=1)
            ax.axvline(start - 0.5, color='red', linestyle='-', alpha=0.8, linewidth=1)
    
    # Add non-standard residue boundaries in RED
    for res_name, res_idx, start, end in non_standard_ranges:
        ax.axhline(start - 0.5, color='red', linestyle='--', alpha=0.8, linewidth=0.5)
        ax.axhline(end - 0.5, color='red', linestyle='--', alpha=0.8, linewidth=0.5)
        ax.axvline(start - 0.5, color='red', linestyle='--', alpha=0.8, linewidth=0.5)
        ax.axvline(end - 0.5, color='red', linestyle='--', alpha=0.8, linewidth=0.5)
    
    # Set axes properties
    ax.set_xlim(-0.5, total_length - 0.5)
    ax.set_ylim(-0.5, total_length - 0.5)
    ax.set_xlabel("Scored residue", fontsize=11)
    ax.set_ylabel("Aligned residue", fontsize=11)
    
    # Compact label positioning
    label_offset_close = total_length * 0.025   # Non-standard residues (closer)
    label_offset_far = total_length * 0.055     # Chain labels (further)
    
    # Add non-standard residue labels first (closer to plot)
    for res_name, res_idx, start, end in non_standard_ranges:
        mid_pos = (start + end) / 2
        # Top
        ax.text(mid_pos, total_length + label_offset_close, f"{res_name}{res_idx}", 
                ha='center', va='bottom', fontsize=8,
                bbox=dict(boxstyle="round,pad=0.2", facecolor="lightblue", alpha=0.8))
        # Right
        ax.text(total_length + label_offset_close, mid_pos, f"{res_name}{res_idx}", 
                ha='center', va='center', fontsize=8, rotation=-90,
                bbox=dict(boxstyle="round,pad=0.2", facecolor="lightblue", alpha=0.8))
    
    # Add chain labels (further from plot)
    for chain_id, mol_label, start, end in chain_ranges:
        mid_pos = (start + end) / 2
        # Top
        ax.text(mid_pos, total_length + label_offset_far, f"{chain_id}\n{mol_label}", 
                ha='center', va='bottom', fontsize=9,
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8, edgecolor='gray'))
        # Right
        ax.text(total_length + label_offset_far, mid_pos, f"{chain_id}\n{mol_label}", 
                ha='center', va='center', fontsize=9, rotation=-90,
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8, edgecolor='gray'))
    
    # Add horizontal colorbar with compact positioning
    cbar = fig.colorbar(im, ax=ax, orientation='horizontal', 
                       shrink=0.7, pad=0.08, aspect=30)
    cbar.set_label(cbar_label, fontsize=10)
    cbar.ax.tick_params(labelsize=9)
    
    # Move title higher to avoid overlap with labels
    fig.suptitle(f"{matrix_type.upper()} - Model {model_num}", 
                fontsize=14, fontweight='bold', x=0.5, y=0.98, ha='center')
    
    # Multi-line explanatory text at bottom
    fig.text(0.5, 0.01, explanation, ha='center', va='bottom', fontsize=8,
             bbox=dict(boxstyle="round,pad=0.3", facecolor="lightgray", alpha=0.8),
             transform=fig.transFigure)
    
    # Compact layout with minimal margins
    plt.subplots_adjust(left=0.08, bottom=0.15, right=0.85, top=0.93)
    
    # Save with tight bounding box to minimize white space
    plt.savefig(os.path.join(output_dir, f'{matrix_type}_model_{model_num}.png'), 
               dpi=300, bbox_inches='tight', pad_inches=0.1)
    plt.close()


def plot_all_models_plddt(directory, runname, output_dir, chain_ranges, non_standard_ranges, total_length):
    """Plot pLDDT values for all models in one figure"""
    plddt_files = sorted(glob(os.path.join(directory, f'plddt_{runname}_model_*.npz')))
    plt.figure(figsize=(14, 7))
    add_confidence_background()
    for i, file in enumerate(plddt_files):
        data = np.load(file)
        plddt = data['plddt'][:total_length] * 100
        plt.plot(range(len(plddt)), plddt, alpha=0.7, label=f'Model {i}', linewidth=1.5)
    add_chain_labels_1d(chain_ranges, non_standard_ranges)
    add_confidence_legend()
    plt.xlim(0, total_length)
    plt.ylim(0, 100)
    plt.ylabel("pLDDT Confidence")
    plt.xlabel("Residue/Atom number")
    plt.title("pLDDT Confidence Scores - All Models")
    plt.legend(loc='upper left')
    explanation = (
        "pLDDT (predicted Local Distance Difference Test) estimates per-residue confidence.\n"
        "Background bands show confidence ranges: >90 (very high), 70-90 (confident), "
        "50-70 (low), <50 (very low). Higher scores indicate more reliable structure predictions."
    )
    plt.figtext(0.5, 0.01, explanation, ha='center', va='bottom', fontsize=9,
                bbox=dict(boxstyle="round,pad=0.5", facecolor="lightgray", alpha=0.8))
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.15)
    plt.savefig(os.path.join(output_dir, f'{runname}_plddt_all_models.png'), dpi=300, bbox_inches='tight')
    plt.close()

def plot_single_model_plddt(directory, runname, output_dir, chain_ranges, non_standard_ranges, total_length, model_num):
    """Plot pLDDT for a single model with AlphaFold colored line and shading"""
    file = os.path.join(directory, f'plddt_{runname}_model_{model_num}.npz')
    data = np.load(file)
    plddt = data['plddt'][:total_length] * 100
    plt.figure(figsize=(14, 6))
    cmap = create_alphafold_colormap()
    x_vals = np.arange(len(plddt))
    for i in range(len(plddt)):
        color = cmap(plddt[i]/100)
        plt.fill_between([i, i+1], [0, 0], [plddt[i], plddt[i]], color=color, alpha=0.3, linewidth=0)
    points = np.array([x_vals, plddt]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    colors = [cmap(plddt[i]/100) for i in range(len(plddt)-1)]
    lc = LineCollection(segments, colors=colors, linewidths=1.5, alpha=0.8)
    plt.gca().add_collection(lc)
    add_chain_labels_1d(chain_ranges, non_standard_ranges)
    add_confidence_legend()
    plt.xlim(0, total_length)
    plt.ylim(0, 100)
    plt.ylabel("pLDDT Confidence")
    plt.xlabel("Residue/Atom number")
    plt.title(f"pLDDT Confidence Scores - Model {model_num}")
    explanation = (
        "Line and shading colors indicate confidence levels at each position. "
        "Blue regions show high confidence predictions, while red/orange indicate uncertainty.\n"
        "Chain boundaries marked with grey lines, non-standard residues with dashed lines."
    )
    plt.figtext(0.5, 0.01, explanation, ha='center', va='bottom', fontsize=9,
                bbox=dict(boxstyle="round,pad=0.5", facecolor="lightgray", alpha=0.8), wrap=True)
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.15)
    plt.savefig(os.path.join(output_dir, f'{runname}_plddt_model_{model_num}.png'), dpi=300, bbox_inches='tight')
    plt.close()

def process_confidence_data(directory, runname, model_num, chain_id_map):
    """Process confidence JSON file and return metadata"""
    json_file = os.path.join(directory, f'confidence_{runname}_model_{model_num}.json')
    if not os.path.exists(json_file):
        return None, None
    
    with open(json_file) as f:
        data = json.load(f)
    
    # Process metadata for CSV (same as before)
    metadata = {
        'model': model_num,
        'confidence_score': round(data.get('confidence_score', 0), 3),
        'ptm': round(data.get('ptm', 0), 3),
        'iptm': round(data.get('iptm', 0), 3),
        'ligand_iptm': round(data.get('ligand_iptm', 0), 3),
        'complex_plddt': round(data.get('complex_plddt', 0), 3),
        'complex_iplddt': round(data.get('complex_iplddt', 0), 3)
    }
    
    # Process chains_ptm
    chains_ptm = data.get('chains_ptm', {})
    for eid, value in chains_ptm.items():
        if float(value) == 0.0:
            continue
        chain_info = chain_id_map.get(int(eid))
        if chain_info:
            metadata[f'ptm_{chain_info["id"]}'] = round(float(value), 3)
    
    # Create matrix sized by number of chains, not residues
    num_chains = len(chain_id_map)
    iptm_matrix = np.zeros((num_chains, num_chains), dtype=np.float32)
    
    # Process pair_chains_iptm
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
                
            # Add to metadata
            metadata[f'iptm_{src_chain["id"]}_{dst_chain["id"]}'] = round(value, 3)
            
            # Add to matrix using chain indices, not residue indices
            iptm_matrix[src_chain['index'], dst_chain['index']] = value
    
    return metadata, iptm_matrix

def plot_iptm_matrix(directory, runname, output_dir, model_num, chain_ranges, iptm_matrix):
    """Plot chain interaction ipTM matrix with interface-focused visualization"""
    plt.figure(figsize=(8, 7))
    chain_ids = [f"{chain_id}\n{mol_label}" for chain_id, mol_label, _, _ in chain_ranges]
    num_chains = len(chain_ids)
    symmetric_matrix = np.maximum(iptm_matrix, iptm_matrix.T)
    masked_matrix = np.ma.masked_where(symmetric_matrix == 0, symmetric_matrix)
    cmap = plt.get_cmap('viridis')
    im = plt.imshow(masked_matrix, cmap=cmap, origin='lower', vmin=0, vmax=1, aspect='equal')
    plt.xticks(range(num_chains), chain_ids, fontsize=10)
    plt.yticks(range(num_chains), chain_ids, fontsize=10)
    cbar = plt.colorbar(im, label='Interface Confidence (ipTM)', shrink=0.8)
    cbar.set_ticks(np.arange(0, 1.1, 0.2))
    for i in range(num_chains):
        for j in range(num_chains):
            if symmetric_matrix[i, j] > 0:
                plt.text(j, i, f'{symmetric_matrix[i, j]:.2f}', 
                        ha='center', va='center', 
                        color='white' if symmetric_matrix[i, j] < 0.5 else 'black',
                        fontsize=9)
            else:
                plt.text(j, i, 'N/A', ha='center', va='center', color='grey', fontsize=8, style='italic')
    plt.title(f"Protein Interface Confidence - Model {model_num}")
    plt.xlabel("Interacting Chain")
    plt.ylabel("Interacting Chain")
    explanation = (
        "ipTM scores estimate interface prediction confidence (0-1).\n"
        "Values show mutual interface quality between chain pairs.\n"
        "Higher scores indicate more reliable interactions."
        )
    plt.figtext(0.5, 0.01, explanation, ha='center', va='bottom', fontsize=9,
                bbox=dict(boxstyle="round,pad=0.5", facecolor="lightgray", alpha=0.8), wrap=True)
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.18)
    plt.savefig(os.path.join(output_dir, f'{runname}_iptm_model_{model_num}.png'), dpi=300, bbox_inches='tight')
    plt.close()

def main():
    parser = argparse.ArgumentParser(description='Generate Boltz-2 analysis plots')
    parser.add_argument('input_dir', help='Path to directory containing Boltz-2 results')
    parser.add_argument('-o', '--output', default='plots', help='Output directory')
    args = parser.parse_args()
    
    os.makedirs(args.output, exist_ok=True)
    runname = detect_runname(args.input_dir)
    chains, residues = load_chains_data(args.input_dir, runname)
    chain_ranges, non_standard_ranges, total_length = parse_chain_ranges_with_residues(chains, residues)
    
    # Create chain ID mapping
    chain_id_map = {}
    for idx, (chain_id, mol_label, start, end) in enumerate(chain_ranges):
        chain_data = next(c for c in chains if c['name'].item() == chain_id)
        chain_id_map[chain_data['asym_id'].item()] = {
            'id': chain_id,
            'type': mol_label,
            'index': idx
        }
    
    # CSV metadata collection
    all_metadata = []
    
    # Plot summary pLDDT plot (only once, outside model loop)
    plot_all_models_plddt(args.input_dir, runname, args.output, chain_ranges, non_standard_ranges, total_length)
    
    # Find number of models
    num_models = len(glob(os.path.join(args.input_dir, f'plddt_{runname}_model_*.npz')))
    
    for model_num in range(num_models):
        # Per-model pLDDT plot
        plot_single_model_plddt(args.input_dir, runname, args.output, chain_ranges, non_standard_ranges, total_length, model_num)
        
        # PAE and PDE plots (if files exist)
        if os.path.exists(os.path.join(args.input_dir, f'pae_{runname}_model_{model_num}.npz')):
            plot_2d_matrix(args.input_dir, runname, args.output, 'pae', model_num, chain_ranges, non_standard_ranges, total_length)
        if os.path.exists(os.path.join(args.input_dir, f'pde_{runname}_model_{model_num}.npz')):
            plot_2d_matrix(args.input_dir, runname, args.output, 'pde', model_num, chain_ranges, non_standard_ranges, total_length)
        
        # Process confidence data and ipTM plot
        metadata, iptm_matrix = process_confidence_data(args.input_dir, runname, model_num, chain_id_map)
        if metadata and iptm_matrix is not None:
            all_metadata.append(metadata)
            plot_iptm_matrix(args.input_dir, runname, args.output, model_num, chain_ranges, iptm_matrix)
    
    # Write CSV metadata
    if all_metadata:
        csv_path = os.path.join(args.output, f'{runname}_metadata.csv')
        # Collect all unique fieldnames
        fieldnames = set()
        for md in all_metadata:
            fieldnames.update(md.keys())

        with open(csv_path, 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(all_metadata)
    
    # Write separate affinity CSV (if present)
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