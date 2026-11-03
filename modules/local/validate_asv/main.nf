process VALIDATE_ASV {

    tag "ASV input validation"

    container 'asv-explorer:0.2.0'

    publishDir "${params.outdir}", mode: 'copy'

    input:
    path counts
    path taxonomy

    output:
    path "validation.json"
    path "sample_summary.tsv"
    path "feature_summary.tsv"

    script:
    """
    python /app/scripts/validate_asv.py \
        --counts "${counts}" \
        --taxonomy "${taxonomy}" \
        --outdir .
    """
}