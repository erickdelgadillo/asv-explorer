process VALIDATE_ASV {

    tag "ASV input validation"

    publishDir "${params.outdir}", mode: 'copy'

    input:
    path counts
    path taxonomy

    output:
    path "validation_preview.txt"

    script:
    """
    echo "Counts file: ${counts.name}" > validation_preview.txt
    echo "Taxonomy file: ${taxonomy.name}" >> validation_preview.txt
    """
}