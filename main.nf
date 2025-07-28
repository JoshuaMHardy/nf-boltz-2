#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    main:
    
    if(params.input_yaml) {
        input_yaml_path = file(params.input_yaml)
        
        if(input_yaml_path.isDirectory()) {
            println "Using existing YAMLs from: ${params.input_yaml}"
            
            // Create channel from YAML files in directory
            yaml_files = Channel.fromPath("${params.input_yaml}/*.yaml")
                .map { yaml_file -> 
                    def pred_name = yaml_file.getSimpleName()
                    return [yaml_file, pred_name]
                }
            
            yaml_files.view { "Found YAML: ${it[1]}" }
            
        } else if(input_yaml_path.isFile()) {
            println "Using single YAML file: ${params.input_yaml}"
            pred_name = input_yaml_path.getSimpleName()
            
            // Create channel with single file
            yaml_files = Channel.of([input_yaml_path, pred_name])
            
        } else {
            error "Input YAML path does not exist: ${params.input_yaml}"
        }
        
    } else {
        println "Generating YAML from parameters"
        
        // Generate YAML from parameters
        GenerateYAML(params)
        
        // Create channel from generated YAML
        yaml_files = GenerateYAML.out.yaml.map { yaml_file -> 
            return [yaml_file, params.pred_name]
        }
    }
    
    // Run Boltz for each YAML file
    RunBoltz(yaml_files)

    // Generate metadata and plots
    AnalyseBoltz(RunBoltz.out.boltz_results)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PROCESSES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process RunBoltz {
    container "/vast/projects/alphafold/databases/predictedDb/containers/proteinDJ/21_05/boltz.sif"
    containerOptions "--nv --bind ${params.boltz_models}:/boltzcache"

    queue 'gpuq'
    clusterOptions '--gres=gpu:' + params.gpu_model + ':1'
    cpus params.n_cpus
    memory params.memory
    executor 'slurm'

    publishDir "${params.outdir}", mode: 'copy', pattern: "${pred_name}"
    publishDir "${params.outdir}", mode: 'copy', pattern: "*.log"

    input:
    tuple path(yaml), val(pred_name)

    output:
    path "${pred_name}"
    tuple path("boltz_results_*/processed/structures/*.npz"), path("boltz_results_*/predictions/*/*"), val(pred_name), emit: boltz_results
    path ("*.log")

    script:
    """
    # We specify a tmp directory as some python packages try to write to the user home directory outside the container
    mkdir tmp
    export HOME=tmp

    boltz predict \
        ${yaml} \
        --cache /boltzcache \
        --use_msa_server \
        ${params.diffusion_samples ? "--diffusion_samples ${params.diffusion_samples}" : ''} \
        ${params.method ? "--method ${params.method}" : ''} \
        ${params.recycling_steps ? "--recycling_steps ${params.recycling_steps}" : ''} \
        ${params.sampling_steps ? "--sampling_steps ${params.sampling_steps}" : ''} \
        ${params.affinity_mw_correction ? '--affinity_mw_correction ' : ''} \
        ${params.diffusion_samples_affinity ? "--diffusion_samples_affinity ${params.diffusion_samples_affinity}" : ''} \
        ${params.sampling_steps_affinity ? "--sampling_steps ${params.sampling_steps_affinity}" : ''} \
        ${params.output_format ? "--output_format ${params.output_format} " :'' } \
        ${params.write_full_pae ? '--write_full_pae ' : ''} \
        ${params.write_full_pde ? '--write_full_pde ' : ''} \
        ${params.extra_config ? params.boltz_extra_config : ''} \
        2>&1 | tee boltz_${pred_name}.log

    # Create prediction directory
    mkdir ${pred_name}
   
    # Copy predictions and input yaml to prediction directory
    cp boltz_results_*/predictions/*/*.${params.boltz_output_format} ${pred_name}/.
    cp ${yaml} ${pred_name}/.
    """
}

process AnalyseBoltz {
    container "/vast/projects/alphafold/databases/predictedDb/containers/proteinDJ/21_05/pyrosetta_tools.sif"
    containerOptions "-B ${projectDir}/bin/generate_metadata.py:/generate_metadata.py"

    cpus 4
    memory '8G'
    executor 'slurm'

    publishDir "${params.outdir}/${pred_name}/", mode: 'copy'

    input:
    tuple path(struct_npz_file), path(model_npz_files), val(pred_name)

    output:
    path "*.png"
    path "*.csv"

    script:
    """
    eval "\$(micromamba shell hook --shell bash)"
    micromamba activate pyrosetta

    # Generate plots and metadata for Boltz-1 prediction 
    python /generate_metadata.py "./" -o "./"

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
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
