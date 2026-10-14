library(reticulate)
use_condaenv("~/anaconda3/envs/r4.3/")
library(Seurat)
library(dplyr)
library(gplots)
library(argparse)
library(cellxgene.census)
library(data.tree)
set.seed(123)


parser <- ArgumentParser(description = "Process Seurat objects and transfer labels.")

parser$add_argument("--integration_method", type="character", help="Integration method of query and reference", default="pcaproject")
parser$add_argument("--ref_keys", type="character", nargs="*", help="List of reference keys to pass to query_transfer", default=c("subclass", "class", "family"))
parser$add_argument("--dims", type="integer", help="Number of dimensions", default=50)
parser$add_argument("--max.features", type="integer", help="Maximum number of features", default=200)
parser$add_argument("--k.anchor", type="integer", help="Number of anchors", default=10)
parser$add_argument("--k.score", type="integer", help="?", default=30)
parser$add_argument("--cutoff", type="numeric", help="Cutoff threshold for label transfer prediction scores", default=0)
parser$add_argument("--ref_path", type="character", help="path to references")
parser$add_argument("--query_path", type="character", help="path to query")
parser$add_argument("--k.weight", type="integer", help="k.weight", default=50)
parser$add_argument("--normalization_method", type="character", help="Normalization method", default="SCT")
parser$add_argument("--nfeatures", type="integer", help="Number of variable features to use for dim reduction", default=2000)


args <- parser$parse_args()
if (is.null(args$ref_path) || is.null(args$query_path)) stop("--ref_path and --query_path are required")
# Extract arguments from the parsed list    
#batch_key <- args$batch_key
integration_method <- args$integration_method
ref_keys <- args$ref_keys
dims <- args$dims
max.features <- args$max.features
k.anchor <- args$k.anchor
k.score <- args$k.score
project.query <- NULL
k.weight <- args$k.weight
ref_path <- args$ref_path
query_path <- args$query_path
normalization_method <- args$normalization_method
n_features <- args$nfeatures

ref = readRDS(ref_path)
query = readRDS(query_path)

transfer_multiple_labels <- function(
            query, reference, reduction, normalization_method="SCT",
            ref_keys, dims, max.features, 
            k.anchor, k.score, k.weight, project.query=NULL) {
    
    #threshold_dims <- 2000
    if (is.null(project.query)) {
        if ((ncol(query) > ncol(reference))) {
            project.query = TRUE
        } else {
            project.query = FALSE
        }
    }


    anchors <- FindTransferAnchors(reference=reference, 
        normalization.method=normalization_method, 
        query=query, 
        # assay should default to whatever normalization was applied
        npcs=dims, dims=1:dims, 
        reduction = reduction, # had to remove project query for SCTransform to work
        max.features=max.features, k.anchor=k.anchor, k.score=k.score,
        reference.reduction="pca") # use precomputed PCA from seurat_processing step


# need to fix this 
    k.weight = min(k.weight, floor(nrow(anchors@anchors) / k.score ))

    key = ref_keys[1] # assumes keys are ordered

    predictions <- tryCatch({
    
        TransferData(anchorset = anchors, refdata=key, reference=reference, weight.reduction=reduction, k.weight = k.weight)

        }, error = function(e) {
        
        message("Error in TransferData: ", e$message)
        # Extract k.weight suggestion (e.g., "less than 27") from error message
        # this only works if the error message format is consistent
        # check if the error message contains "less than"
        if (grepl("less than", e$message)) {
            match <- regmatches(e$message, regexec("less than ([0-9]+)", e$message))
            k.weight.new <- as.integer(match[[1]][2]) - 1 # Decrease by 1 to ensure it is less than the suggested value
            message("Adjusting k.weight to: ", k.weight.new)
            # Retry TransferData with the corrected k.weight
            TransferData(
                anchorset = anchors,
                refdata = key,
                reference = reference,
                weight.reduction = reduction,
                k.weight = k.weight.new
            )
        } else {
            message("Unable to adjust k.weight from error message. Setting k.weight to FALSE.")
            # set k.weight to FALSE
            TransferData(
                anchorset = anchors,
                refdata = key,
                reference = reference,
                weight.reduction = reduction,
                k.weight = FALSE
            )
        }
    })
    return(predictions)

}


prediction_scores <- transfer_multiple_labels(normalization_method=normalization_method,
        query=query, reference=ref, reduction=integration_method, 
        ref_keys=ref_keys, dims=dims, 
        max.features=max.features, 
        k.anchor=k.anchor, k.score=k.score, 
        project.query=project.query, k.weight=k.weight)



query_name = basename(query_path) %>% gsub("_processed.rds", "", .)
ref_name = basename(ref_path) %>% gsub(".rds", "", .)

# Seurat rewrites class names in the score columns (make.names style). Map each column back to the reference label by
# matching on lower-case alphanumerics only, and stop if any column or label cannot be matched uniquely.
norm_label <- function(x) tolower(gsub("[^A-Za-z0-9]", "", x))
ref_labels <- unique(as.character(ref@meta.data[[ref_keys[1]]]))
prediction_scores <- prediction_scores %>%
  as.data.frame() %>%
  select(-any_of(c("predicted.id", "prediction.score.max")))
score_labels <- sub("^prediction\\.score\\.", "", colnames(prediction_scores))
idx <- match(norm_label(score_labels), norm_label(ref_labels))
if (anyNA(idx) || anyDuplicated(norm_label(ref_labels))) {
  stop("Could not map Seurat score columns to reference labels: ", paste(score_labels[is.na(idx)], collapse = ", "))
}
colnames(prediction_scores) <- ref_labels[idx]

write.table(prediction_scores, file=gzfile(paste0(query_name,"_",ref_name,"_prediction_scores_seurat.tsv.gz"), compression=4), sep="\t", row.names=FALSE)
write.table(query@meta.data, file=paste0(query_name,".obs.relabel.tsv"), row.names=FALSE, sep= "\t")