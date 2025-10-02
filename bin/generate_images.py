#!/usr/bin/env python

import argparse
import sys
from chimerax.core.commands import run

# Parse command line arguments
parser = argparse.ArgumentParser(description='Generate ChimeraX images for protein structures')
parser.add_argument('-binder', '--binder', type=str, help='Binder chain ID')
parser.add_argument('--extension', type=str, default='pdb', help='File extension (pdb or cif)')

args = parser.parse_args()
binder_id = args.binder
extension = args.extension

# Open structure files
run(session, f'open *.{extension}')

# Set up graphics
run(session, 'window 512 512')
run(session, 'preset "overall look" "publication 2 (depth-cued)"')
run(session, 'graphics selection color black width 2')

# Align to first model
for model in session.models.list():
    id = int(model.id_string)
    if id > 1:
        run(session, f'matchmaker #{id} to #1')

# Hide models for now
run(session, 'hide models')

for model in session.models.list():
    # Add color key (model ID is created but not in the list above)
    run(session, 'color bfactor palette alphafold key true')
    run(session, 'key fontSize 9 size 0.5,0.04000 pos 0.25,0.05')
    
    # Remove file extension from name
    model_name = model.name.split('.')[0]
    id = int(model.id_string)

    run(session, f'show #{id} models')
    
    # Select and highlight binder if present
    if binder_id:
        run(session, f'select #{id}/{binder_id}')

    # Front image
    run(session, f'view #{id}')
    run(session, f'2dlabels text "{model_name}" xpos 0.3 ypos 0.95 size 12 color black')
    run(session, f'save "{model_name}_overview_front.png" width 1024 height 1024')
    run(session, '2dlabels delete all')

    # Side image
    run(session, 'turn y 90')
    run(session, f'view #{id}')
    run(session, f'2dlabels text "{model_name}" xpos 0.3 ypos 0.95 size 12 color black')
    run(session, f'save "{model_name}_overview_side.png" width 1024 height 1024')
    run(session, '2dlabels delete all')
    run(session, 'turn y -90')

    # Binder zoom images (only if binder_id is provided)
    if binder_id:
        print(f"Processing binder chain: {binder_id}")
        
        # Select binder chain
        run(session, f'select #{id}/{binder_id}')
        
        # Make everything else semi-transparent
        run(session, f'select ~#{id}/{binder_id} & #{id}')
        run(session, f'transparency sel 70 atoms,cartoons')
        
        # Clear selection and select binder again
        run(session, 'select clear')
        run(session, f'select #{id}/{binder_id}')
        
        # Zoom to binder
        run(session, f'view #{id}/{binder_id}')
        
        # Find hydrogen bonds and label/show residues
        run(session, f'hbonds #{id}/{binder_id} restrict #{id} intraMol false intraRes false select true reveal true')
        run(session, f'label sel text "{{0.chain_id}} {{0.name}} {{0.number}}{{0.insertion_code}}"')
        
        # Front binder image
        run(session, f'2dlabels text "{model_name}" xpos 0.3 ypos 0.95 size 12 color black')
        run(session, f'select #{id}/{binder_id}')  # Highlight binder
        run(session, f'save "{model_name}_binder_front.png" width 1024 height 1024')
        run(session, '2dlabels delete all')

        # Side binder image
        run(session, 'turn y 90')
        run(session, f'2dlabels text "{model_name}" xpos 0.3 ypos 0.95 size 12 color black')
        run(session, f'save "{model_name}_binder_side.png" width 1024 height 1024')
        run(session, '2dlabels delete all')
        run(session, 'turn y -90')
        
        # Reset transparency and clear selections
        run(session, f'transparency #{id} 0 atoms,cartoons')
        run(session, 'select clear')
        run(session, 'label delete')
    else:
        print("No binder ID provided, skipping binder images")

    # Hide model for next round of images
    run(session, f'hide #{id} models')

run(session, 'exit')
