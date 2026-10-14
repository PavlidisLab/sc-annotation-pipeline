process REF_PROCESS_SEURAT {
    tag "$h5ad_file.baseName"
    label 'process_high'

    conda "/home/rschwartz/anaconda3/envs/r4.3"

    input:
    path h5ad_file

    output:
    path "${h5ad_file.baseName}.rds", emit: ref_rds

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    Rscript ${projectDir}/bin/ref_preprocessing.R \\
        --h5ad_file ${h5ad_file} \\
        --normalization_method ${params.seurat_normalization} \\
        --dims ${params.seurat_dims} \\
        --nfeatures ${params.seurat_nfeatures} \\
        --k.score ${params.seurat_k_score} \\
        --k.anchor ${params.seurat_k_anchor} \\
        --k.weight ${params.seurat_k_weight} \\
        ${args}
    """

    stub:
    """
    touch ${h5ad_file.baseName}.rds
    """
}
