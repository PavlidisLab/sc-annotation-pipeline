process PREDICT_SEURAT {
    tag "$query_name"
    label 'process_high'

    conda "/home/rschwartz/anaconda3/envs/r4.3"

    input:
    tuple val(study_name), val(query_name), path(query_h5ad), path(query_rds), path(ref_rds)
    val ref_keys

    output:
    tuple val(study_name), val(query_name), path(query_h5ad), path("*_prediction_scores_seurat.tsv.gz"), emit: scores

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    Rscript ${projectDir}/bin/predict_seurat.R \\
        --query_path ${query_rds} \\
        --ref_path ${ref_rds} \\
        --ref_keys ${ref_keys[0]} \\
        --integration_method pcaproject \\
        --dims ${params.seurat_dims} \\
        --max.features ${params.seurat_max_features} \\
        --k.score ${params.seurat_k_score} \\
        --k.anchor ${params.seurat_k_anchor} \\
        --k.weight ${params.seurat_k_weight} \\
        --normalization_method ${params.seurat_normalization} \\
        ${args}
    """

    stub:
    """
    touch ${query_name}_ref_prediction_scores_seurat.tsv.gz
    """
}
