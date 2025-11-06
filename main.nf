#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    validateInputs()

    // Calculate output directory once - outdir can be relative or absolute path
    def output_dir = params.use_outdir_prefix ? joinPaths(params.outdir_prefix, params.outdir) : params.outdir
    full_outdir = file(output_dir)

    // Log the output directory configuration
    if (params.use_outdir_prefix) {
        println("Outdir_prefix: '${params.outdir_prefix}'")
    }
    println("Output directory: '${output_dir}'")

    if (params.input_yaml) {
        input_yaml_path = file(params.input_yaml)

        if (input_yaml_path.isDirectory()) {
            println("Using existing YAMLs from: ${params.input_yaml}")

            // Create channel from YAML files in directory
            yaml_files = channel
                .fromPath("${params.input_yaml}/*.yaml")
                .map { yaml_file ->
                    def pred_name = yaml_file.getSimpleName()
                    return [yaml_file, pred_name]
                }

            yaml_files.view { "Found YAML: ${it[1]}" }

            template_files_ch = []
        }
        else if (input_yaml_path.isFile()) {
            println("Using single YAML file: ${params.input_yaml}")
            pred_name = input_yaml_path.getSimpleName()

            // Create channel with single file
            yaml_files = channel.of([input_yaml_path, pred_name])

            template_files_ch = []
        }
        else {
            error("Input YAML path does not exist: ${params.input_yaml}")
        }
    }
    else {
        // Check for templates
        def (template_files_list, pdb_codes_list) = collectTemplateFilesAndPDBs(params)

        files_ch = channel.fromList(template_files_list)
        cif_files_ch = channel.empty()

        if (pdb_codes_list.size() > 0) {
            println("Downloading templates")
            DownloadTemplates(pdb_codes_list)
            cif_files_ch = DownloadTemplates.out.cif_files
        }

        // Merge all files to a single channel for staging
        template_files_ch = files_ch.mix(cif_files_ch).collect().ifEmpty([])

        println("Generating YAML from parameters")
        GenerateYAML(params)
        yaml_files = GenerateYAML.out.yaml.map { yaml_file ->
            [yaml_file, params.pred_name]
        }
    }
    // Check if we need to run ColabFold for MSA
    PrepColabFold(yaml_files)
    
    // Collect any pre-existing MSA files referenced in YAMLs
    CollectPreexistingMSAs(yaml_files)

    // Combine CSV files from PrepColabFold, keeping only one header
    PrepColabFold.out
        .collectFile(
            name: 'combined_msa.csv',
            keepHeader: true,
            skip: 1,
        )
        .set { seqs_for_msa }
    PrepColabFold.out.ifEmpty { log.info("Skipping ColabFold - no MSAs to generate") }

    // Run ColabFold only if there are CSV files
    RunColabFold(seqs_for_msa, output_dir)

    // Combine generated MSAs from ColabFold with pre-existing MSAs from YAML
    msa_files_ch = RunColabFold.out.csv
        .mix(CollectPreexistingMSAs.out.msa_files.flatten())
        .collect()
        .ifEmpty([])

    // Run Boltz for each YAML file
    RunBoltz(yaml_files, template_files_ch, msa_files_ch, output_dir)

    // Generate metadata and plots
    AnalyseBoltz(RunBoltz.out.boltz_results, output_dir)

    // Generate ChimeraX plots
    VisualiseBoltz(RunBoltz.out.boltz_models, output_dir)

    // Generate PDF of plots
    //GeneratePDF(AnalyseBoltz.out.plots,VisualiseBoltz.out)

    // Save log file on completion
    def out_pred_name = params.pred_name

    workflow.onComplete {
        def logFile = file('.nextflow.log')
        if (logFile.exists()) {
            logFile.copyTo(full_outdir.resolve("${out_pred_name}/nextflow.log"))
        }
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PROCESSES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
process GenerateYAML {
    input:
    val params

    output:
    path "*.yaml", emit: yaml

    script:
    """
    generate_yaml.py '${groovy.json.JsonOutput.toJson(params)}' > ${params.pred_name}.yaml
    """
}

process DownloadTemplates {
    input:
    val pdb_codes_tuples

    output:
    path "*.cif", emit: cif_files

    when:
    pdb_codes_tuples.size() > 0

    script:
    def download_cmds = pdb_codes_tuples
        .collect { pdb_code, filename ->
            "wget -q -O ${filename} https://files.rcsb.org/download/${pdb_code}.cif"
        }
        .join('\n')
    """
    ${download_cmds}
    for file in *.cif; do
      if [ ! -s "\$file" ]; then
        echo "Error: Could not download \$file"
        exit 1
      fi
    done
    """
}

process CollectPreexistingMSAs {
    // Scan YAML files for pre-existing MSA file paths
    input:
    tuple path(yaml), val(pred_name)

    output:
    path ("msa_files/*"), optional: true, emit: msa_files

    script:
    """
    #!/usr/bin/env python3
    import yaml
    import os
    import shutil
    
    # Read and parse the YAML file
    with open('${yaml}', 'r') as f:
        data = yaml.safe_load(f)
    
    # Create output directory
    os.makedirs('msa_files', exist_ok=True)
    
    # Collect any MSA file paths from the YAML
    for block in data.get('sequences', []):
        if 'protein' in block:
            protein = block['protein']
            msa = protein.get('msa')
            
            # Check if msa is a file path (not True, False, or "empty")
            if msa and isinstance(msa, str) and msa not in ['empty', 'true', 'false']:
                # This is a file path - check if it exists and copy it
                if os.path.exists(msa):
                    basename = os.path.basename(msa)
                    shutil.copy(msa, f'msa_files/{basename}')
                    print(f"Found pre-existing MSA: {msa}")
                else:
                    print(f"Warning: MSA file specified but not found: {msa}")
    """
}

process PrepColabFold {
    // Check YAML files to see if we need to generate MSA(s)
    input:
    tuple path(yaml), val(pred_name)

    output:
    path ("*.csv"), optional: true, emit: csv

    script:
    """
    #!/usr/bin/env python3
    import yaml
    
    # Read and parse the YAML file
    with open('${yaml}', 'r') as f:
        data = yaml.safe_load(f)
    
    sequences = []
    
    # Process each protein block
    for block in data.get('sequences', []):
        if 'protein' in block:
            protein = block['protein']
            sequence = protein.get('sequence')
            msa = protein.get('msa')
            
            # Check if we need to generate MSA for this sequence
            # Skip if: msa is "empty", or msa is a PRE-EXISTING file path
            # Include if: msa is a future filename (no path separator) or True
            needs_generation = False
            
            if msa == 'empty':
                # Explicitly no MSA
                needs_generation = False
            elif isinstance(msa, str) and '/' in msa:
                # Has path separator - this is a pre-existing file path, skip generation
                needs_generation = False
            elif msa:
                # msa is True or a filename we need to generate (no path separator)
                needs_generation = True
            
            if needs_generation and sequence:
                # Add one copy per chain ID
                chain_ids = protein.get('id', [])
                if isinstance(chain_ids, str):
                    chain_ids = [chain_ids]
                for chain_id in chain_ids:
                    sequences.append(sequence)
    
    # Write CSV if we have sequences
    if sequences:
        seq_string = ':'.join(sequences)
        csv_content = f"id,sequence\\n${pred_name},{seq_string}\\n"
        
        with open('${pred_name}.csv', 'w') as f:
            f.write(csv_content)
    """
}



process RunColabFold {
    container "${params.colabfold_container}"
    containerOptions "--nv --bind ${params.colabfold_dbs}:/colabfold_dbs --bind ${projectDir}/bin/split_msa.py:/split_msa.py"

    cpus params.colabfold_cpus
    memory params.colabfold_memory
    queue 'gpuq'
    clusterOptions '--gres=gpu:' + params.colabfold_gpu_model + ':1'

    publishDir "${output_dir}/msas", mode: 'copy'

    input:
    path csv
    val output_dir

    output:
    path ("split_msas/*.csv"), emit: csv

    script:

    """
    # Run ColabFold on combined CSV file to generate MSA per row
    colabfold_search \
        --db-load-mode 2 \
        --gpu 1 \
        --prefilter-mode 1 \
        --threads ${task.cpus} \
        ${csv} \
        /colabfold_dbs \
        "results/"

    # Split each MSA by unique query sequence for Boltz-2 input
    python /split_msa.py -i results/ -o split_msas
    """
}

process RunBoltz {
    container "${params.boltz2_container}"
    containerOptions "--nv --bind ${params.boltz2_models}:/boltzcache"

    queue 'gpuq'
    clusterOptions '--gres=gpu:' + params.boltz2_gpu_model + ':1'
    cpus params.boltz2_cpus
    memory params.boltz2_memory

    publishDir "${output_dir}", mode: 'copy', pattern: "${pred_name}"

    input:
    tuple path(yaml), val(pred_name)
    path template_files
    path msa_files
    val output_dir

    output:
    path "${pred_name}"
    tuple path("boltz_results_*/processed/structures/*.npz"), path("boltz_results_*/predictions/*/*"), val(pred_name), emit: boltz_results
    tuple path("${pred_name}/*"), val(pred_name), emit: boltz_models

    script:
    """
    # We specify a tmp directory as some python packages try to write to the user home directory outside the container
    mkdir tmp
    export HOME=tmp
    
    # Create prediction directory to copy results into
    mkdir ${pred_name}

    boltz predict \
        ${yaml} \
        --cache /boltzcache \
        ${params.diffusion_samples ? "--diffusion_samples ${params.diffusion_samples}" : ''} \
        ${params.method ? "--method '${params.method}'" : ''} \
        ${params.use_potentials ? '--use_potentials' : ''} \
        ${params.recycling_steps ? "--recycling_steps ${params.recycling_steps}" : ''} \
        ${params.sampling_steps ? "--sampling_steps ${params.sampling_steps}" : ''} \
        ${params.step_scale ? "--step_scale ${params.step_scale}" : ''} \
        ${params.affinity_mw_correction ? '--affinity_mw_correction ' : ''} \
        ${params.diffusion_samples_affinity ? "--diffusion_samples_affinity ${params.diffusion_samples_affinity}" : ''} \
        ${params.sampling_steps_affinity ? "--sampling_steps ${params.sampling_steps_affinity}" : ''} \
        ${params.output_format == "pdb" ? "--output_format pdb " : ''} \
        ${params.write_full_pae ? '--write_full_pae ' : ''} \
        ${params.write_full_pde ? '--write_full_pde ' : ''} \
        ${params.chunk_size_transition_z ? "--chunk_size_transition_z ${params.chunk_size_transition_z}" : ''} \
        ${params.chunk_size_transition_msa ? "--chunk_size_transition_msa ${params.chunk_size_transition_msa}" : ''} \
        ${params.chunk_size_tri_attn ? "--chunk_size_tri_attn ${params.chunk_size_tri_attn}" : ''} \
        ${params.triangle_mult_gate_nchunks ? "--triangle_mult_gate_nchunks ${params.triangle_mult_gate_nchunks}" : ''} \
        ${params.chunk_size_outer_product ? "--chunk_size_outer_product ${params.chunk_size_outer_product}" : ''} \
        ${params.chunk_size_threshold ? "--chunk_size_threshold ${params.chunk_size_threshold}" : ''} \
        ${params.use_bfloat16 ? '--use_bfloat16' : ''} \
        ${params.extra_config ? "'${params.extra_config}'" : ''} \
        2>&1 | tee ${pred_name}/boltz.log
   
    # Copy predictions and input yaml to prediction directory. Use .cif if mmcif format.
    cp boltz_results_*/predictions/*/*.${params.output_format == 'mmcif' ? 'cif' : 'pdb'} ${pred_name}/.
    cp ${yaml} ${pred_name}/.
    # Copy MSA files if they exist
    for file in *.csv; do
        [ -e "\$file" ] || continue
        cp "\$file" ${pred_name}/.
    done

    """
}

process AnalyseBoltz {
    container "/vast/projects/alphafold/databases/predictedDb/containers/proteinDJ/21_05/pyrosetta_tools.sif"
    containerOptions "-B ${projectDir}/bin/generate_metadata.py:/generate_metadata.py"

    cpus 4
    memory '8G'

    publishDir "${output_dir}/${pred_name}/", mode: 'copy', pattern: "*.csv"
    publishDir "${output_dir}/${pred_name}/plots", mode: 'copy', pattern: "*.png"

    input:
    tuple path(struct_npz_file), path(model_npz_files), val(pred_name)
    val output_dir

    output:
    path "*.png", emit: plots
    path "*.csv"

    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate pyrosetta

    # Generate plots and metadata for Boltz-2 prediction 
    python /generate_metadata.py "./" -o "./"

    """
}

process VisualiseBoltz {
    container "/stornext/System/data/software/rhel/9/base/structbio/chimerax/1.10.1/chimerax-1.10.1.sif"
    containerOptions "--writable-tmpfs -B /home -B ${projectDir}/bin/generate_images.py:/generate_images.py"

    cpus 4
    memory '24G'

    publishDir "${output_dir}/${pred_name}/images", mode: 'copy', pattern: "*.png"

    input:
    tuple path(pred_files), val(pred_name)
    val output_dir

    output:
    tuple path("*.png"), val(pred_name)

    script:
    """
    # Extract file extension for script
    extension="${params.output_format == 'mmcif' ? 'cif' : 'pdb'}"
    
    # Check for binder in YAML and extract chain ID
    binder_arg=""
    if grep -q "affinity:" *.yaml && grep -q "binder:" *.yaml; then
        binder=\$(awk '/binder:/ {gsub(/^[[:space:]]*/, "", \$2); gsub(/[[:space:]]*\$/, "", \$2); print \$2}' *.yaml)
        if [[ -n "\$binder" ]]; then
            echo "Found binder chain ID: \$binder"
            binder_arg="--binder \$binder"
        fi
    else
        echo "No binder found in YAML file"
    fi

    # Run ChimeraX with appropriate arguments
    chimerax --nogui --offscreen --script "/generate_images.py --extension \$extension \$binder_arg"
    
    # List generated images
    ls -la *.png
    """
}


process GeneratePDF {
    cpus 4
    memory '8G'
    maxRetries 2

    publishDir "${output_dir}/${pred_name}/", mode: 'copy', pattern: "*.pdf"

    input:
    path plots
    tuple path(images), val(pred_name)
    val output_dir

    output:
    path "*.pdf"

    script:
    """
    module load ImageMagick
    model_numbers=\$(ls *_model_*.png | sed 's/.*_model_\\([0-9]\\+\\)_.*/\\1/' | sort -u)
    for model in \$model_numbers; do
        echo "Processing model \$model..."
        magick *_\${model}_plddt.png *_\${model}_pae.png *_\${model}_pde.png *_\${model}_iptm.png -background white ${pred_name}_model_\${model}_plots.pdf
        
        # Check if binder images exist
        binder_images=\$(ls *_\${model}_binder*.png 2>/dev/null)
        if [ -n "\$binder_images" ] ; then
            montage *_\${model}_overview*.png \$binder_images -tile 2x2 -geometry 512x512+32+32 -border 2 -background white ${pred_name}_model_\${model}_images.pdf
        else
            montage *_\${model}_overview*.png -tile 2x1 -geometry 512x512+32+32 -border 2 -background white ${pred_name}_model_\${model}_images.pdf
        fi
    done
    """
}



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
def validateInputs() {
    if (params.input_yaml && !file(params.input_yaml).exists()) {
        exit(1, "ERROR: Input YAML file not found: ${params.input_yaml}")
    }

    if (!params.seq1_seq && !params.input_yaml) {
        exit(1, "ERROR: At least one sequence must be provided")
    }

    // Validate prediction name
    if (!params.pred_name.matches(/^[a-zA-Z0-9_-]+$/)) {
        exit(1, "ERROR: pred_name can only contain letters, numbers, underscores, and dashes")
    }
}

def joinPaths(prefix, suffix) {
    // Smart path joining that handles slashes correctly
    if (!prefix || prefix.isEmpty()) {
        return suffix
    }
    if (!suffix || suffix.isEmpty()) {
        return prefix
    }

    // Remove trailing slash from prefix if present
    def cleanPrefix = prefix.endsWith('/') ? prefix[0..-2] : prefix
    // Remove leading slash from suffix if present
    def cleanSuffix = suffix.startsWith('/') ? suffix[1..-1] : suffix

    // Join with single slash
    return "${cleanPrefix}/${cleanSuffix}"
}

def collectTemplateFilesAndPDBs(params) {
    def files = []
    def pdb_codes = []

    // Handle explicit template1/2 fields
    [1, 2].each { i ->
        def path = params["template${i}_path"]
        def pdb_code = params["template${i}_pdb_code"]
        if (pdb_code) {
            pdb_codes << [pdb_code.toUpperCase(), "${pdb_code.toUpperCase()}.cif"]
        }
        else if (path) {
            files << file(path)
        }
    }

    // Handle template_extra (assume only PDB codes, no paths)
    if (params.template_extra) {
        // regex: [XXXX] with XXXX = 4-char alphanum code
        def re = /\[([0-9A-Za-z]{4})\]/
        (params.template_extra =~ re).each { match ->
            def code = match[1].toUpperCase()
            pdb_codes << [code, "${code}.cif"]
        }
    }

    return [files.unique(), pdb_codes.unique()]
}
