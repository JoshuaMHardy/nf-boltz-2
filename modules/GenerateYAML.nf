process GenerateYAML {
    input:
    val params
    
    output:
    path 'input.yaml', emit: yaml
    
    script:
    """
    generate_yaml.py '${groovy.json.JsonOutput.toJson(params)}' > input.yaml
    """
}