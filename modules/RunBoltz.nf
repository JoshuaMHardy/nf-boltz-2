process RunBoltz {
    container "/vast/projects/alphafold/databases/predictedDb/containers/proteinDJ/21_05/boltz.sif"
    containerOptions "--nv"

    queue 'gpuq'
    clusterOptions '--gres=gpu:' + params.gpu_model + ':1'
    cpus params.n_cpus
    memory params.memory
    executor 'slurm'

    publishDir "${params.out_dir}", mode: 'copy'

    input:
    path yaml

    output:
    path "${params.pred_name}"
    path ("*.log"), emit: logs

    script:
    """
        # We specify tmp directories as some python packages try to write to the user home directory outside the container
        mkdir tmp
        export NUMBA_CACHE_DIR=tmp
        export XDG_CONFIG_HOME=tmp
        export TRITON_CACHE_DIR=tmp

        # Rename input yaml to prediction name
        mv input.yaml ${params.pred_name}.yaml

        boltz predict \
            ${params.pred_name}.yaml \
            --output_format pdb \
            --diffusion_samples ${params.boltz_diffusion_samples} \
            --recycling_steps ${params.boltz_recycling_steps} \
            --sampling_steps ${params.boltz_sampling_steps} \
            --cache /boltz/cache \
            --use_msa_server \
            ${params.boltz_extra_config ? params.boltz_extra_config : ''} \
            2>&1 | tee ${params.pred_name}.log

        # Move prediction directory
        mv boltz_results_*/predictions/${params.pred_name} .
        # Copy YAML to prediction directory
        cp ${params.pred_name}.yaml ${params.pred_name}/.
        """
}
