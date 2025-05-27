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