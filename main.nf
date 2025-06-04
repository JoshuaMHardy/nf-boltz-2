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
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PROCESSES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process RunBoltz {
    container "/vast/projects/alphafold/databases/predictedDb/containers/proteinDJ/21_05/boltz.sif"
    containerOptions "--nv"

    queue 'gpuq'
    clusterOptions '--gres=gpu:' + params.gpu_model + ':1'
    cpus params.n_cpus
    memory params.memory
    executor 'slurm'

    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple path(yaml), val(pred_name)

    output:
    path "${pred_name}"
    path ("*.log"), emit: logs

    script:
    """
    # We specify tmp directories as some python packages try to write to the user home directory outside the container
    mkdir tmp
    export NUMBA_CACHE_DIR=tmp
    export XDG_CONFIG_HOME=tmp
    export TRITON_CACHE_DIR=tmp

    boltz predict \
        ${yaml} \
        --output_format pdb \
        --diffusion_samples ${params.boltz_diffusion_samples} \
        --recycling_steps ${params.boltz_recycling_steps} \
        --sampling_steps ${params.boltz_sampling_steps} \
        --cache /boltz/cache \
        --use_msa_server \
        ${params.boltz_write_full_pae ? '--write_full_pae ' : ''} \
        ${params.boltz_write_full_pde ? '--write_full_pde ' : ''} \
        ${params.boltz_extra_config ? params.boltz_extra_config : ''} \
        2>&1 | tee boltz_${pred_name}.log

    # Move prediction directory
    mv boltz_results_*/predictions/${pred_name} .
    # Copy input yaml to prediction directory
    cp ${yaml} ${pred_name}/.
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
