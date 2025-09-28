#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    validateInputs()

    if (params.input_yaml) {
        input_yaml_path = file(params.input_yaml)

        if (input_yaml_path.isDirectory()) {
            println("Using existing YAMLs from: ${params.input_yaml}")

            // Create channel from YAML files in directory
            yaml_files = Channel
                .fromPath("${params.input_yaml}/*.yaml")
                .map { yaml_file ->
                    def pred_name = yaml_file.getSimpleName()
                    return [yaml_file, pred_name]
                }

            yaml_files.view { "Found YAML: ${it[1]}" }
        }
        else if (input_yaml_path.isFile()) {
            println("Using single YAML file: ${params.input_yaml}")
            pred_name = input_yaml_path.getSimpleName()

            // Create channel with single file
            yaml_files = Channel.of([input_yaml_path, pred_name])
        }
        else {
            error("Input YAML path does not exist: ${params.input_yaml}")
        }
    }
    else {
        // Check for templates
        def (template_files_list, pdb_codes_list) = collectTemplateFilesAndPDBs(params)

        files_ch = Channel.fromList(template_files_list)
        cif_files_ch = Channel.empty()

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
    // Run Boltz for each YAML file
    RunBoltz(yaml_files, template_files_ch)

    // Generate metadata and plots
    AnalyseBoltz(RunBoltz.out.boltz_results)

    // Generate ChimeraX plots
    VisualiseBoltz(RunBoltz.out.boltz_models)

    // Generate PDF of plots
    GeneratePDF(AnalyseBoltz.out.plots,VisualiseBoltz.out)

    // Save log file on completion
    def out_pred_name = params.pred_name
    def outputDir = file(params.outdir)

    workflow.onComplete {
        def logFile = file('.nextflow.log')
        if (logFile.exists()) {
            logFile.copyTo(outputDir.resolve("${out_pred_name}/nextflow.log"))
        }
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PROCESSES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
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

process RunBoltz {
    container "${params.container_path}"
    containerOptions "--nv --bind ${params.models_path}:/boltzcache"

    queue 'gpuq'
    clusterOptions '--gres=gpu:' + params.gpu_model + ':1'
    cpus params.n_cpus
    memory params.memory
    executor 'slurm'

    publishDir "${params.outdir}", mode: 'copy', pattern: "${pred_name}"

    input:
    tuple path(yaml), val(pred_name)
    path template_files

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
        --use_msa_server \
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
        ${params.extra_config ? "'${params.extra_config}'" : ''} \
        2>&1 | tee ${pred_name}/boltz.log
   
    # Copy predictions and input yaml to prediction directory. Use .cif if mmcif format.
    cp boltz_results_*/predictions/*/*.${params.output_format == 'mmcif' ? 'cif' : 'pdb'} ${pred_name}/.
    cp ${yaml} ${pred_name}/.
    """
}

process AnalyseBoltz {
    container "/vast/projects/alphafold/databases/predictedDb/containers/proteinDJ/21_05/pyrosetta_tools.sif"
    containerOptions "-B ${projectDir}/bin/generate_metadata.py:/generate_metadata.py"

    cpus 4
    memory '8G'
    executor 'slurm'

    publishDir "${params.outdir}/${pred_name}/", mode: 'copy', pattern: "*.csv"
    publishDir "${params.outdir}/${pred_name}/plots", mode: 'copy', pattern: "*.png"

    input:
    tuple path(struct_npz_file), path(model_npz_files), val(pred_name)

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
    cpus 4
    memory '8G'
    executor 'slurm'

    publishDir "${params.outdir}/${pred_name}/images", mode: 'copy', pattern: "*.png"

    input:
    tuple path(pred_files), val(pred_name)

    output:
    tuple path("*.png"), val(pred_name)

    script:
    """
    # Prepare files to remove soft-links incompatible with ChimeraX
    cp ${projectDir}/bin/generate_images.py .
    mkdir inputs
    cp -f *.${params.output_format == 'mmcif' ? 'cif' : 'pdb'} inputs/.

    module load chimerax
    
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
    chimerax --nogui --offscreen --script "generate_images.py --extension \$extension \$binder_arg"
    
    # List generated images
    ls -la *.png
    """
}


process GeneratePDF {
    cpus 4
    memory '8G'
    executor 'slurm'
    
    publishDir "${params.outdir}/${pred_name}/", mode: 'copy', pattern: "*.pdf"

    input:
    path(plots)
    tuple path(images), val(pred_name)

    output:
    path "*.pdf"

    script:
    """
    module load ImageMagick
    model_numbers=\$(ls *_model_*.png | sed 's/.*_model_\\([0-9]\\+\\)_.*/\\1/' | sort -u)
    for model in \$model_numbers; do
        echo "Processing model \$model..."
        magick *_\${model}_plddt.png *_\${model}_pae.png *_\${model}_pde.png *_\${model}_iptm.png -background white ${pred_name}_model_\${model}_plots.pdf
        montage *_\${model}_overview*.png *_\${model}_binder*.png -tile 2x2 -geometry 512x512+32+32 -border 2 -background white ${pred_name}_model_\${model}_images.pdf
    done
    """
}

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
