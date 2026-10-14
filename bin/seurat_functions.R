load_required_packages <- function() {
  packages <- c("Seurat", "dplyr", "ggplot2", "data.table", "readxl", "stringr", "harmony","leiden", 
  "gprofiler2", "GEOquery", "ggExtra","SeuratDisk","tidyr","patchwork","RColorBrewer")

  for (pkg in packages) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
  
}

process_sample <- function(raw_counts, sample) {
    tryCatch({
        seurat_obj <- CreateSeuratObject(counts = raw_counts, project = sample, min.cells = 3, min.features = 200)
        seurat_obj@meta.data$sample <- sample
        seurat_obj@meta.data$orig.ident <- rownames(seurat_obj@meta.data)
        seurat_obj <- SetIdent(seurat_obj, value= seurat_obj$sample)
       # seurat_obj_list[[sample]] <- seurat_obj
    }, error = function(e) {
        print(paste("All cells in sample ", sample, "fail filtering thresholds, moving on to next"))
        print(e)
     }
    )
    return(seurat_obj)
}

load_seurat <- function(topdir, mode) {
    seurat_obj_list <- list()

    if (mode %in% c("directory", "files")) {
        # Initialize an empty list to store Seurat objects
        
        
        if (mode == "directory") {
            subdirs <- list.dirs(topdir, full.names = TRUE, recursive = FALSE)
            samples <- sapply(strsplit(subdirs, "/"), function(x) { #Splits each path in subdirs by /
            #Extracts the last component of each split result and splits by .
            #extracts first component
              parts <- strsplit(tail(x, 1), "\\.")[[1]] 
              return(parts)
            })

        } else {  # mode == "files"
            filenames <- dir(topdir, pattern = "(barcodes\\.tsv|features\\.tsv|\\.mtx)")

            samples <- unique(sub("_([^_]+)$", "", filenames))
            #samples <- unique(sub("barcodes\\.*|features\\.*|\\.mtx\\.*", "", filenames))
            #samples <- unique(sapply(strsplit(filenames, "_"), function(x) { 
               #  return(x)[[1]]
            #}))
        }

        for (sample in samples) {
            if (mode == "directory") {
                raw_counts <- Read10X(file.path(topdir, sample))
            } else {  # mode == "files"
                matching_files <- dir(topdir, pattern = sample, full.names = TRUE)
                barcodes_file <- matching_files[grep('barcodes', matching_files)]
                features_file <- matching_files[grep('features|genes', matching_files)]
                matrix_file <- matching_files[grep('mtx', matching_files)]
                raw_counts <- ReadMtx(mtx = matrix_file, features = features_file, cells = barcodes_file)
            }
            seurat_obj_list[[sample]] = process_sample(raw_counts, sample)
            
        }
        
        print(names(seurat_obj_list))
        seurat_obj <- merge(seurat_obj_list[[1]], seurat_obj_list[2:length(seurat_obj_list)], add.cell.ids =names(seurat_obj_list))#project = "seurat_obj")
        seurat_obj <- JoinLayers(seurat_obj)
        
        return(seurat_obj)
    } else {
        stop("Invalid mode. Mode should be 'directory' or 'files'")
    }
}

# adapted from https://matthieuxmoreau.github.io/EarlyPallialNeurogenesis/html-Reports/Quality_Control.html

plot_prefilter <- function(outdir, dataset) {
    # Visualize QC metrics as a violin plot
    options(repr.plot.width = 30, repr.plot.height = 20)
    plot0 <- VlnPlot(dataset, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), ncol = 4, pt.size = 0.01)
    plot1 <- FeatureScatter(dataset, feature1 = "nCount_RNA", feature2 = "percent.mt")
    plot2 <- FeatureScatter(dataset, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

    # Combine plots using cowplot
    combined_plot <- plot0 / plot1 / plot2
    return(combined_plot)
    #ggsave(filename=paste0(outdir,"_prefilter.png"),width=20,height=20)

}
  

get_mt_threshold <- function(dataset, nmads) {
    Cell.QC.Stat <- dataset@meta.data
    max.mito.thr <- median(Cell.QC.Stat$percent.mt) + nmads*mad(Cell.QC.Stat$percent.mt)
    min.mito.thr <- median(Cell.QC.Stat$percent.mt) - nmads*mad(Cell.QC.Stat$percent.mt)
    return(list("max"=max.mito.thr,"min"=min.mito.thr))
}
       
get_rp_threshold <- function(dataset, nmads) {
    Cell.QC.Stat <- dataset@meta.data
    max.rp.thr <- median(Cell.QC.Stat$percent.rp) + nmads*mad(Cell.QC.Stat$percent.rp)
    min.rp.thr <- median(Cell.QC.Stat$percent.rp) - nmads*mad(Cell.QC.Stat$percent.rp)
    return(list("max"=max.rp.thr,"min"=min.rp.thr))
}

plot_mito_rp <- function(outdir, dataset, nmads, base_theme) {
    Cell.QC.Stat <- dataset@meta.data
    #Filtering cells based on percentage of mitochondrial transcripts
    #We applied a high and low median absolute deviation (mad) thresholds to exclude outlier cells

    max.mito.thr <- get_mt_threshold(dataset, nmads$mito)[["max"]]
    min.mito.thr <- get_mt_threshold(dataset, nmads$mito)[["min"]]

    Cell.QC.Stat$valid_mt<-  with(Cell.QC.Stat, 
    percent.mt <= max.mito.thr & 
    percent.mt >= min.mito.thr
    )

    mt <- ggplot(Cell.QC.Stat, aes(x=nFeature_RNA, y=percent.mt, color=valid_mt)) +
          geom_point() +
          geom_hline(aes(yintercept = max.mito.thr), colour = "red", linetype = 2) +
          geom_hline(aes(yintercept = min.mito.thr), colour = "red", linetype = 2) +
          base_theme #+
         # annotate(geom = "text", 
                   # label = paste0(as.numeric(table(Cell.QC.Stat$valid_mt)[2])," cells removed"), 
                  #  x = 6000, y = 0)

    mtplot <- ggMarginal(mt, type = "histogram", fill="lightgrey", bins=100, size=2) 

    max.rp.thr <- get_rp_threshold(dataset, nmads$rp)[["max"]]
    min.rp.thr <- get_rp_threshold(dataset, nmads$rp)[["min"]]

    Cell.QC.Stat$valid_rp<-  with(Cell.QC.Stat, 
        percent.rp <= max.rp.thr & 
        percent.rp >= min.rp.thr
        )

    rp <- ggplot(Cell.QC.Stat, aes(x=nFeature_RNA, y=percent.rp, color=valid_rp)) +
          geom_point() +
          base_theme +
          geom_hline(aes(yintercept = max.rp.thr), colour = "red", linetype = 2) +
          geom_hline(aes(yintercept = min.rp.thr), colour = "red", linetype = 2) 
          #annotate(geom = "text", label = paste0(as.numeric(table(Cell.QC.Stat$valid_rp)[2])," cells removed"), x = 6000, y = 0)


    rpplot <- ggMarginal(rp, type = "histogram", fill="lightgrey", bins=100, size=2)

    mt_rp <- patchwork::wrap_plots(mtplot, rpplot, nrow = 1) & base_theme
    return(mt_rp)
   # ggsave(filename=paste0(outdir, "mt_rp.png"), plot=mt_rp)
}

get_genes_threshold <- function(dataset, nmads) {
    Cell.QC.Stat <- dataset@meta.data
    # Set low and hight thresholds on the number of detected genes
    min.Genes.thr <- median(log10(Cell.QC.Stat$nFeature_RNA)) - nmads*mad(log10(Cell.QC.Stat$nFeature_RNA))
    max.Genes.thr <- median(log10(Cell.QC.Stat$nFeature_RNA)) + nmads*mad(log10(Cell.QC.Stat$nFeature_RNA))
    return(list("max"=max.Genes.thr,"min"=min.Genes.thr)) 
    }

get_umi_threshold <-  function(dataset, nmads) {
    Cell.QC.Stat <- dataset@meta.data
    # Set max threshold on the number of transcripts
    max.nUMI.thr <- median(log10(Cell.QC.Stat$nCount_RNA)) + nmads*mad(log10(Cell.QC.Stat$nCount_RNA))
    min.nUMI.thr <- median(log10(Cell.QC.Stat$nCount_RNA)) - nmads*mad(log10(Cell.QC.Stat$nCount_RNA))
    return(list("max"=max.nUMI.thr,"min"=min.nUMI.thr)) 
    }

get_lm <- function(dataset, nmads=4) {
    Cell.QC.Stat <- dataset@meta.data
    # Gene/UMI scatter plot before filtering
    lm.model <- lm(data = Cell.QC.Stat, formula = log10(nFeature_RNA) ~ log10(nCount_RNA))

    # Calculate residuals
    residuals <- resid(lm.model)
    # Calculate MAD of the residuals
    mad_residuals <- mad(residuals)
    # Calculate the intercept adjustment using MAD
    intercept_adjustment <- quantile(residuals, 0.5) - nmads * mad_residuals  # Adjust MAD factor based on your preference
    return(list("model"=lm.model,"intercept_adjustment"=intercept_adjustment))
}

plot_lmfit_pre <- function(outdir, dataset, nmads, base_theme) {
    
    Cell.QC.Stat <- dataset@meta.data

    # Set low and hight thresholds on the number of detected genes
    min.Genes.thr <- get_genes_threshold(dataset,nmads$genes)[["min"]]
    max.Genes.thr <- get_genes_threshold(dataset,nmads$genes)[["max"]]

    # Set max threshold on the number of transcripts
    max.nUMI.thr <- get_umi_threshold(dataset,nmads$umi)[["max"]]
    min.nUMI.thr <- get_umi_threshold(dataset,nmads$umi)[["min"]]

    lm <- get_lm(dataset, nmads$lm)
    lm.model <- lm[["model"]]
    intercept_adjustment <- lm[["intercept_adjustment"]]


    Cell.QC.Stat$valid_cells <-  with(Cell.QC.Stat, 
    log10(nFeature_RNA) > min.Genes.thr & 
    log10(nCount_RNA) < max.nUMI.thr & 
    log10(nCount_RNA) > min.nUMI.thr & 
    log10(nFeature_RNA) < max.Genes.thr &
    log10(nFeature_RNA) > (log10(nCount_RNA) * lm.model$coefficients[2] + (lm.model$coefficients[1] + intercept_adjustment)) &
    log10(nFeature_RNA) < (log10(nCount_RNA) * lm.model$coefficients[2] + (lm.model$coefficients[1] - intercept_adjustment)))
    
    plot5 <- ggplot(Cell.QC.Stat, aes(x=log10(nCount_RNA), y=log10(nFeature_RNA), color=valid_cells)) +
      geom_point() +
      base_theme +
      geom_smooth(method="lm") +
      geom_hline(aes(yintercept = min.Genes.thr), colour = "green", linetype = 2) +
      geom_hline(aes(yintercept = max.Genes.thr), colour = "green", linetype = 2) +
      geom_vline(aes(xintercept = max.nUMI.thr), colour = "red", linetype = 2) +
      geom_vline(aes(xintercept = min.nUMI.thr), colour = "red", linetype = 2) +
      geom_abline(intercept = lm.model$coefficients[1] + intercept_adjustment , slope = lm.model$coefficients[2], color="orange")
      #annotate(geom = "text", label = paste0(as.numeric(table(Cell.QC.Stat$valid_cells)[1])," cells removed"), x = 6000, y = 0)
 

    lmfit_pre <- ggMarginal(plot5, type = "histogram", fill="lightgrey")
    return(lmfit_pre)
}

get_valid_cells <- function(dataset, nmads=list(mito=4, rp=4, genes=4, umi=4, lm=4)) {
    Cell.QC.Stat <- dataset@meta.data
    # Set low and hight thresholds on the number of detected genes
    max.mito.thr <- get_mt_threshold(dataset,nmads$mito)[["max"]]
    min.mito.thr <- get_mt_threshold(dataset,nmads$mito)[["min"]]
    
    max.rp.thr <- get_rp_threshold(dataset,nmads$rp)[["max"]]
    min.rp.thr <- get_rp_threshold(dataset,nmads$rp)[["min"]]
    
    min.Genes.thr <- get_genes_threshold(dataset,nmads$genes)[["min"]]
    max.Genes.thr <- get_genes_threshold(dataset,nmads$genes)[["max"]]

    # Set max threshold on the number of transcripts
    max.nUMI.thr <- get_umi_threshold(dataset,nmads$umi)[["max"]]
    min.nUMI.thr <- get_umi_threshold(dataset,nmads$umi)[["min"]]
    
    lm <- get_lm(dataset, nmads$lm)
    lm.model <- lm[["model"]]
    intercept_adjustment <- lm[["intercept_adjustment"]]
     
    # Mark cells as TRUE or FALSE based on the QC thresholds
    Cell.QC.Stat$passes_QC <- with(Cell.QC.Stat, 
    percent.mt <= max.mito.thr & 
    percent.mt >= min.mito.thr & 
    percent.rp <= max.rp.thr & 
    percent.rp >= min.rp.thr & 
    log10(nFeature_RNA) > min.Genes.thr & 
    log10(nCount_RNA) < max.nUMI.thr & 
    log10(nCount_RNA) > min.nUMI.thr & 
    log10(nFeature_RNA) < max.Genes.thr &
    log10(nFeature_RNA) > (log10(nCount_RNA) * lm.model$coefficients[2] + (lm.model$coefficients[1] + intercept_adjustment)) &
    log10(nFeature_RNA) < (log10(nCount_RNA) * lm.model$coefficients[2] + (lm.model$coefficients[1] - intercept_adjustment)
    ))

    # Return the modified data frame
    return(Cell.QC.Stat)

}


filter_valid_cells <- function(dataset) {
    Cell.QC.Stat <- dataset@meta.data
    #Cell.QC.Stat <- get_valid_cells(dataset, nmads=nmads)
    Cell.QC.Stat <- Cell.QC.Stat %>% filter(Cell.QC.Stat$passes_QC)
    filtered <- dataset[ ,rownames(Cell.QC.Stat)]
    return(filtered)
}

plot_post_filtering <- function(outdir, dataset) {

    plot7 <- VlnPlot(dataset, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), ncol = 4, pt.size = 0.01)
    plot8 <- FeatureScatter(dataset, feature1 = "nCount_RNA", feature2 = "percent.mt")
    plot9 <- FeatureScatter(dataset, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    combined_plot <- plot7 / plot8 / plot9
    return(combined_plot)
    #ggsave(filename=paste0(outdir,"_postfilter.png"), plot=combined_plot, width=20,height=20)
    
}

plot_qc_metrics <- function(outdir, dataset, nmads=list(mito=4, rp=4, genes=4, umi=4, lm=4)) {
# adapted from https://matthieuxmoreau.github.io/EarlyPallialNeurogenesis/html-Reports/Quality_Control.html
    base_theme <- theme(
        text = element_text(size = 25),            # Base text size
        axis.title = element_text(size = 25),      # Axis title size
        axis.text = element_text(size = 25),       # Axis tick label size
        plot.title = element_text(size = 25, face = "bold")  # Plot title size
    )
    # Generate plots
    p2 <- plot_mito_rp(outdir, dataset, nmads=nmads, base_theme)
    p3 <- plot_lmfit_pre(outdir, dataset, nmads=nmads, base_theme)

    # Check if malat1_threshold is in dataset@meta.data
    if ("malat1_threshold" %in% colnames(dataset@meta.data)) {
        p4 <- DimPlot(dataset, group.by = c("malat1_threshold", "passes_QC")) & base_theme

    } else {
        p4 <- DimPlot(dataset, group.by = "passes_QC") & base_theme
    }
    # Combine plots using patchwork
    combined_plot <- p2 /p3/ p4 + base_theme

    # Save the combined plot
    ggsave(file.path(outdir, "qc_metrics_patchwork.png"), combined_plot,height=30, width=20)
}

load_cell_cycle_genes <- function(organism) {
    if (organism=="mus_musculus") {
       s.genes = gorth(cc.genes.updated.2019$s.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
       g2m.genes = gorth(cc.genes.updated.2019$g2m.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
    }
    
    if (organism=="homo_sapiens") {
        s.genes <- cc.genes$s.genes
        g2m.genes <- cc.genes$g2m.genes
    }
    geneslist <- list()
    geneslist[["s.genes"]] <- s.genes
    geneslist[["g2m.genes"]] <- g2m.genes
    
    return(geneslist)
}

score_cell_cycle <- function(seurat_obj, split.key="sample", organism, scale_factor=10000, normalization_method="LogNormalize") {
    
    cell.cycle.genes <- load_cell_cycle_genes(organism=organism)
    
    seurat_obj_list <- SplitObject(seurat_obj, split.by=split.key)
    seurat_obj_list <- lapply(X = seurat_obj_list, FUN = function(x) {
    x <- NormalizeData(x, scale.factor=scale_factor, normalization.method=normalization_method)
    x <- CellCycleScoring(x, s.features = cell.cycle.genes[["s.genes"]], g2m.features = cell.cycle.genes[["g2m.genes"]])
    })
    seurat_obj <- merge(seurat_obj_list[[1]], seurat_obj_list[2:length(seurat_obj_list)])  
    return(seurat_obj)
}

harmony_integration <- function(seurat_obj, nFeatures=2000, dims=30, resolution=0.8, vars.to.regress=NULL, integration_key="sample",
                               scale_factor=10000, normalization_method="LogNormalize", features=NULL) {   
    
    integration_key_value <- eval(parse(text = paste0("seurat_obj$", integration_key)))
    seurat_obj[["RNA"]] <- split(seurat_obj[["RNA"]], f = integration_key_value)

    seurat_obj <- seurat_obj %>% ScaleData(vars.to.regress=vars.to.regress) %>% RunPCA(npcs=dims, reduction_name="pca") 
    seurat_obj <- IntegrateLayers(seurat_obj, method=HarmonyIntegration, orig.reduction = "pca", new.reduction = 'harmony')
    seurat_obj <- JoinLayers(seurat_obj)
    seurat_obj <- seurat_obj %>% RunUMAP(dims=1:dims, reduction="harmony", reduction.name="umap.harmony")
    return(seurat_obj)
}


plot_label_transfer_stats <- function(dataset, outdir) {

    if (!dir.exists(outdir)) {
      # If the directory doesn't exist, create it
      dir.create(dir_path)
      print("Directory created successfully.")
    }    

    density <- ggplot(dataset@meta.data,  aes(x = prediction.score.max, fill = predicted.id)) +
        geom_density(alpha = 0.5) +
        labs(title = "Density Plot of Max Prediction Scores by Cell Type",
        x = "Max Prediction Score", y = "Density") + scale_fill_discrete(name = "Cell Type") +
        facet_wrap(~predicted.id, scales = "free") 

    ggsave(filename=paste0(outdir,"prediction_scores.png", plot=density, width = 20, height = 12))

}

plot_label_stats <- function(dataset, label, sample_key) {
    # Extract identity and sample information from seurat object to determine the number of cells per cluster per sample
    n_cells <- FetchData(dataset, vars = label) %>%
        dplyr::count(!!sym(label))
    colnames(n_cells) <- c(label, "n")

    # Barplot of number of cells per cluster by sample
    p0 <- ggplot(n_cells, aes_string(x = label, y = "n", fill = label)) +
        geom_bar(stat = "identity", position = position_dodge()) +
        geom_text(aes(label = n), vjust = -0.2, position = position_dodge(1)) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) + # Rotate x-axis ticks by 90 degrees 
        labs(title = "Number of Cells per Label", x = label, y = "Number of Cells") +
        guides(fill = guide_legend(title = label)) # Adjust legend title

    # Sample proportions within label
    p1 <- ggplot(dataset@meta.data, aes_string(x = label, fill = sample_key)) +
          geom_bar(position = "fill") +
          theme(axis.text.x = element_text(angle = 90, hjust = 1)) + # Rotate x-axis labels by 90 degrees
          labs(title = "Sample Proportions Within Label", y = "Proportion")

    # Combine plots
    p2 <- p0 + p1
    return(p2)
}


load_geo_metadata <- function(f) {
    acc <- getGEO(f,GSEMatrix=FALSE)
    GSMs <- GSMList(acc)
    meta <- list()
    for (name in names(GSMs)) {
        meta[[name]] <- list()
        gsmlist <- (Meta(GSMs[[name]]))
        chars <- gsmlist$characteristics_ch1
        taxon <- gsmlist$taxid_ch1
        source <- gsmlist$description
        for (char in chars) {
            meta[[name]][[paste0(strsplit(char, ": ")[[1]][1])]] = strsplit(char, ": ")[[1]][2]
        meta[[name]][["taxon"]] = taxon
        meta[[name]][["source"]] = source
    }
    }
    meta_df <- as.data.frame(do.call(rbind, meta))
    meta_df$sample <- rownames(meta_df)
    
    # Convert all columns to factors
    meta_df[] <- lapply(meta_df, function(x) {
        if (is.character(x)) {
            as.factor(x)
        } else {
            x
        }
    })
    
    return(meta_df)
}

add_geo_metadata <- function(GSM, seurat_obj, sample_key="sample") {
    meta_df <- load_geo_metadata(GSM)
    meta <- left_join(seurat_obj@meta.data, meta_df, by=c(sample_key="sample"))
    rownames(meta) <- rownames(seurat_obj@meta.data)
    seurat_obj@meta.data <- meta
    lapply(colnames(seurat_obj@meta.data), function(x) {
        if (is.list(x)) {
            seurat_obj@meta.data[[x]] <- unlist(seurat_obj@meta.data[[x]])
        }
        seurat_obj@meta.data[[x]] <- as.factor(seurat_obj@meta.data[[x]]) 
    })
}

# Function to grep gene patterns
mt_rb_hb_pcs <- function(seurat_obj, dims = 30, pval_cutoff = 0.01) {
    if (!("pca") %in% names(seurat_obj@reductions)) {
        seurat_obj <- seurat_obj %>% RunPCA(npcs=dims)
    }
  # Perform JackStraw analysis if not already done
  if (!("jackstraw" %in% names(seurat_obj@misc))) {
    seurat_obj <- JackStraw(seurat_obj, dims = dims)
  }
  
  # Get significant contaminant genes which contribute to PCs
  sig <- PCASigGenes(seurat_obj, pcs.use = 1:pc_cutoff, pval.cut = pval_cutoff, use.full = FALSE)
  
  # Grep for specific gene patterns
  ribosomal_genes <- grep("^Rpl|^Rps|^Mrp|^RPL|^RPS|^MRP-", sig, value = TRUE)
  mitochondrial_genes <- grep("^mt-|^MT-", sig, value = TRUE)
  hemoglobin_genes <- grep("^Hb|^HB", sig, value = TRUE)
  
  # Return a named list of gene names
  return(c(ribosomal_genes, mitochondrial_genes, hemoglobin_genes))
}

get_mt_rb_percentage <- function(seurat_obj) {

    ribo_genes <- "^Rpl|^Rps|^Mrp|^RPL|^RPS|^MRP"
    mito_genes <- "^mt-|^MT-"
    seurat_obj$percent.mt <- PercentageFeatureSet(seurat_obj, pattern=mito_genes)
    seurat_obj$percent.rp <- PercentageFeatureSet(seurat_obj, pattern=ribo_genes)
    return(seurat_obj)
}

rename_features <- function(seurat_obj, column_name) {
    counts <- seurat_obj[["RNA"]]@counts
    data <- seurat_obj[["RNA"]]@data 

    feature_meta <- seurat_obj[["RNA"]][[]]
    feature_meta[["orig_features"]] <- rownames(feature_meta)
    rownames(feature_meta) <- feature_meta[[column_name]]

    rownames(counts) <- rownames(feature_meta)
    rownames(data) <- rownames(feature_meta)

    newRNA <- CreateAssayObject(counts=counts)
    newRNA$data <- data
    newRNA[[]] <- feature_meta

   seurat_obj[["RNA"]] <- newRNA
   DefaultAssay(seurat_obj) <- "RNA"
   return(seurat_obj)
}

rds_to_h5ad <- function(obj, filepath_prefix) {
    SaveH5Seurat(obj, filename = file.path(filepath_prefix,"h5Seurat"))
    Convert(file.path(filepath_prefix,"h5Seurat"), dest = "h5ad")
}

h5ad_to_rds <- function(h5ad_file_path){
    SeuratDisk::Convert(h5ad_file_path, dest = "h5seurat", overwrite = FALSE)
    message("Reading H5Seurat...")
    h5seurat_file_path <- gsub(".h5ad", ".h5seurat", h5ad_file_path)
    seurat_obj <- SeuratDisk::LoadH5Seurat(h5seurat_file_path, assays = "RNA")
    message("Read Seurat object:")
    return(seurat_obj)
}

malat1_threshold <- function(seurat_obj) {

    malat1_gene <- grep("^(?i)MALAT1$", rownames(seurat_obj[["RNA"]]$data), value = TRUE)
    norm_counts <- seurat_obj[["RNA"]]$data[malat1_gene, ]
    threshold <- define_malat1_threshold(norm_counts) 
    malat1_threshold <- norm_counts > threshold
    seurat_obj$malat1_threshold <- malat1_threshold
    seurat_obj$malat1_threshold <- factor(seurat_obj$malat1_threshold, levels = c("TRUE","FALSE"))
    #good_seurat_obj <- seurat_obj[ ,seurat_obj$malat1_threshold==TRUE]
    return(seurat_obj)
}

process_seurat <- function(seurat_obj, GSE, join_key="sample", vars.to.regress=NULL, organism, outdir, remove_contaminants=FALSE, 
                            filter_mads=TRUE, nmads=list(mito=4, rp=4, genes=4, umi=4, lm=4), filter_malat1=TRUE, dims, score_phase=FALSE, batch_correct=batch_correct, ...) {
    
    if (!is.null(GSE)) {
        # Load GEO metadata and update Seurat object metadata
        meta <- load_geo_metadata(GSE)
        meta <- left_join(seurat_obj@meta.data, meta, by=join_key)
        rownames(meta) <- rownames(seurat_obj@meta.data)
        seurat_obj@meta.data <- meta
    }
    seurat_obj <- get_mt_rb_percentage(seurat_obj)
    if (score_phase) {
        seurat_obj <- score_cell_cycle(seurat_obj, organism=organism)
        seurat_obj$CC.difference <- seurat_obj$G2M.Score - seurat_obj$S.Score
        print("finished scoring cell cycle")
        seurat_obj <- seurat_obj %>% ScaleData() %>% FindVariableFeatures() %>% RunPCA(npcs=dims) %>% RunUMAP(dims=1:dims)
    } else {
        seurat_obj <- seurat_obj %>% NormalizeData() %>% ScaleData() %>% FindVariableFeatures() %>% RunPCA(npcs=dims) %>% RunUMAP(dims=1:dims)
    }
    
    seurat_obj <- tryCatch(
        {
            message("trying to model MALAT1 expression")
            malat1_threshold(seurat_obj)
        },
        error = function(e) {
            message("MALAT1 could not be scored: ", e$message)
            filter_malat1 <<- FALSE  # Use <<- to modify the outer-scope variable
            seurat_obj  # Return the original seurat_obj
        }
    )

    valid_cells <- get_valid_cells(seurat_obj, nmads=nmads)
    seurat_obj@meta.data <- valid_cells 
    plot_qc_metrics(outdir=outdir, dataset=seurat_obj, nmads=nmads)

    if (filter_malat1) { 
       seurat_obj <- seurat_obj[ ,seurat_obj$malat1_threshold==TRUE]     
    }

    if (filter_mads) {
        seurat_obj <- filter_valid_cells(seurat_obj)
    }
    
    if (remove_contaminants) {
        contaminants <- mt_rb_hb_pcs(seurat_obj)
        features=rownames(seurat_obj)[!rownames(seurat_obj) %in% contaminants]
    } else {
        features=NULL
    }
    if (batch_correct) {
        seurat_obj <- harmony_integration(seurat_obj, nFeatures=2000, dims=dims, resolution=0.8, vars.to.regress=vars.to.regress, integration_key=join_key,
                               scale_factor=10000, normalization_method="LogNormalize", features=features)
    }

    return(seurat_obj)
}


process_reference <- function(path_to_reference, scale.factor=10000, normalization.method="LogNormalize", nfeatures=2000) {
    ref <- readRDS(path_to_reference)
    tryCatch({
        counts <- ref[["RNA"]]$counts
    }, warning=function(w) {
        counts <- ref[["RNA"]]$data
    })
    tryCatch({    
        meta <- ref[["RNA"]][[]]
        rownames(counts) <- ref[["RNA"]][[]]$feature_name
        meta$ensembl <- rownames(meta)
        rownames(meta) <- meta$feature_name
        ref_counts <- CreateAssayObject(counts)
        ref[["RNA"]] <- ref_counts
        ref[["RNA"]][[]] <- meta
        DefaultAssay(ref) <- "RNA"
    }, error = function(e) {
        message("no feature_name column, assuming this is not a cellxgene dataset")
    })
    ref <- ref %>% NormalizeData(normalization.method=normalization.method, scale.factor=scale.factor) %>% FindVariableFeatures(nfeatures=nfeatures) %>% ScaleData()
    return(ref)
}

transfer_multiple_labels <- function(
            query, reference, reduction, 
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
    #features = intersect(rownames(reference), rownames(query))
    anchors <- FindTransferAnchors(reference=reference, query=query, npcs=dims, dims=1:dims, reduction = reduction, 
        project.query=project.query, max.features=max.features, k.anchor=k.anchor, k.score=k.score)
    
    k.weight = min(k.weight, floor(nrow(anchors@anchors) / k.score ))
    key = ref_keys[1] # assumes keys are ordered
    #change the k.weight back to 50 or dynamically set ?
    predictions <- TransferData(anchorset = anchors, refdata=key, reference=reference, weight.reduction=reduction, k.weight = k.weight)
    return(predictions)
}


plot_transfer <- function(outdir, reference_name, transfer_method="pcaproject", query, ref_keys, batch_correct=FALSE) {
    if (batch_correct) {
        reduction="umap.harmony"
    } else {
        reduction="umap"
    }

    p <- DimPlot(query, reduction = reduction, group.by = ref_keys, label=FALSE, label.size=6) 
    ggsave(plot=p,filename=paste0(outdir,"/",reference_name,"_label_transfer_",transfer_method,".png"), width=5*length(ref_keys), height=10)
}

get_idents_corr <- function(seurat_obj1, seurat_obj2, assay="RNA") {
  # Get common genes
  common_genes <- intersect(rownames(seurat_obj1[[assay]]$scale.data), rownames(seurat_obj2[[assay]]$scale.data))
  
  # Calculate average expression for common genes
  avg_exp_obj1 <- AggregateExpression(seurat_obj1, assays = assay, features = common_genes)$RNA %>% as.matrix()#, return.seurat=TRUE) %>%
  avg_exp_obj2 <- AggregateExpression(seurat_obj2, assays = assay, features = common_genes)$RNA %>% as.matrix()#, return.seurat=TRUE)
    # Sort the rows by rownames to ensure they are in the same order
  avg_exp_obj1 <- avg_exp_obj1[common_genes, ] 
  avg_exp_obj2 <- avg_exp_obj2[common_genes, ]
  # Convert to numeric matrices
  avg_exp_obj1_cpm <- avg_exp_obj1 %>% cpm(log=TRUE, prior.count=1) 
  avg_exp_obj2_cpm <- avg_exp_obj2 %>% cpm(log=TRUE, prior.count=1)  
  # Compute the correlation matrix
  correlation_matrix <- cor(avg_exp_obj1_cpm, avg_exp_obj2_cpm, use = "pairwise.complete.obs")
  return(correlation_matrix)
}

plot_idents_corr <- function(seurat_obj1, seurat_obj2, assay = "RNA", name1, name2) {
   # Get the names of the datasets

  # Generate the main title using the dataset names
 main_title <- paste("Pearson R", name1, "vs", name2, "")
 
  correlation_matrix <- get_idents_corr(seurat_obj1, seurat_obj2, assay=assay)
  # Plot the heatmap with clustering but without showing the dendrogram
  p <- pheatmap(correlation_matrix,
           #color = colorRampPalette(c("blue", "white", "red"))(100),
           display_numbers = TRUE,
           number_color = "black",
           main = main_title,
           cluster_rows = TRUE,  # Enable clustering of rows
           cluster_cols = TRUE,  # Enable clustering of columns
           show_rownames = TRUE,  # Show row names (optional)
           show_colnames = TRUE,  # Show column names (optional)
           treeheight_row = 0,    # Hide the row dendrogram
           treeheight_col = 0)    # Hide the column dendrogram
  # Convert pheatmap object to ggplot object
  return(as.ggplot(p, scale = 1))
}

highest_corr <- function(seurat_obj1, seurat_obj2, assay="RNA", name1, name2) {
  
  correlation_matrix <- get_idents_corr(seurat_obj1, seurat_obj2, assay=assay) 
  correlation_df <- as.data.frame(as.table(correlation_matrix))
  colnames(correlation_df) <- c(paste(name1), paste(name2), "Correlation")
  #correlation_df$Correlation <- correlation_df$Correlation %>% as.numeric


  # Sort the data frame by the absolute value of correlation in descending order
  highest_correlations <- correlation_df %>%  
    group_by(across(1)) %>% top_n(n=1, Correlation)        
  
  return(highest_correlations)
}

plot_ref_query_corr <- function(query_path, ref_paths, idents, outdir, projPath) { 
  # Load the query dataset
  query <- lapply(query_path, readRDS)
  idents <- lapply(idents, unlist)
  # Load the reference datasets
  refs <- lapply(ref_paths, readRDS)
  #idents <- lapply(idents, unlist)
  # Get the name of the query dataset
  plots <- list()
  query_name <- names(query_path)
  query <- query[[1]]
  # Iterate over each identity for the query dataset
  for (query_ident in idents[[query_name]]) {
    Idents(query) <- eval(parse(text=paste("query", query_ident, sep="$")))
    
    # Iterate over each reference dataset
    for (i in seq_along(refs)) {
      ref_name <- names(refs)[[i]]
      print(ref_name)
      for (ref_ident in idents[[ref_name]]) {
        Idents(refs[[i]]) <- eval(parse(text=paste0("refs[[",i,"]]$", ref_ident)))
        print(ref_ident)
        #print(Idents(refs[[i]]) %>% head())
        # Combine the query with the current reference to create a pair
        pair <- setNames(list(query, refs[[ref_name]]), c(paste(query_name, query_ident, sep=" "), paste(ref_name, ref_ident, sep=" ")))
       # print(pairs)
        
        # Apply the plot_idents_corr function to each pair
       # plots <- lapply(pairs, function(pair) {
          name1 <- names(pair)[[1]]
          name2 <- names(pair)[[2]]
          
          # Create a title by replacing underscores with spaces
          plot_title <- paste0(query_name, "_", query_ident, "_vs_", ref_name, "_", ref_ident)
          
          # Use plot_title in the pheatmap function
          p <- plot_idents_corr(pair[[1]], pair[[2]], name1 = name1, name2 = name2)
          plots[[plot_title]] <- p
          
          x1 <- highest_corr(pair[[1]], pair[[2]], name1 = name1, name2 = name2)
          x2 <- highest_corr(pair[[2]], pair[[1]], name1 = name2, name2 = name1)
          file_name1 <- paste0(query_name, "_", query_ident, "_vs_", ref_name, "_", ref_ident, ".tsv")
          file_name2 <- paste0(ref_name, "_", ref_ident, "_vs_", query_name, "_", query_ident, ".tsv")
          write.table(x1, file = file.path(outdir, file_name1), row.names = FALSE, sep = "\t")
          write.table(x2, file = file.path(outdir, file_name2), row.names = FALSE, sep = "\t") 
        }
      }
    }
  
        

        # Convert pheatmap objects to ggplot objects
        ggplots <- lapply(plots, as.ggplot)
        n_plots <- length(ggplots)
        nrow <- ceiling(sqrt(n_plots))
        ncol <- ceiling(n_plots / nrow)

        # Combine the ggplot objects using patchwork
        grid_plot <- wrap_plots(ggplots, nrow = nrow, ncol = ncol)

        # Save the combined plot
        plot_filename <- paste0(query_name, "_label_compare.png")
        ggsave(file.path(outdir, plot_filename), grid_plot, width = 10 * ncol, height = 5 * nrow)

      }


    relabel <- function(seurat_obj, relabel_df, join_key) {
        meta <- left_join(seurat_obj@meta.data, relabel_df, by=join_key)
        rownames(meta) <- colnames(seurat_obj)
        seurat_obj@meta.data <- meta
        seurat_obj
    }

relabel_wrapper <- function(preproc_obj, relabel_path) {
  # Read relabel table
  relabel_table <- read.table(file = relabel_path, sep = "\t", header = TRUE)
  # Get the join key
  join_key = colnames(relabel_table)[1]
  # Perform relabeling
  relabeled_data <- relabel(preproc_obj, relabel_table, join_key = join_key)
  return(relabeled_data)
}

map_meta <- function(query_transfer, query_file) {
    relabel_path=file.path(projPath,"meta","relabel","human",paste0(query_name,"_relabel.tsv"))
    relabel_table <- read.table(file = relabel_path, sep = "\t", header = TRUE)
    join_key = colnames(relabel_table)[1]
    obj <- query_transfer
    relabel_obj <- relabel(obj, relabel_table, join_key = join_key)
    meta <- relabel_obj@meta.data
    meta$barcode <- colnames(relabel_obj)
    meta
    #write.table(meta, file=outpath, sep="\t", row.names=FALSE)
}

subsample_cells <- function(data, filtered_ids, census, subsample = 500) {
    # Filter data based on filtered_ids
    obs <- data[data$soma_joinid %in% filtered_ids, ]
    celltypes <- unique(obs$cell_type)

    final_idx <- c()

    for (celltype in celltypes) {
        celltype_ids <- obs[obs$cell_type == celltype, ]$soma_joinid
        # Sample if there are enough observations, otherwise take all
        if (length(celltype_ids) >= subsample) {
            subsampled_cell_idx <- sample(celltype_ids, subsample)
        } else {
            subsampled_cell_idx <- celltype_ids
        }
        # Append subsampled indices to final list
        final_idx <- c(final_idx, subsampled_cell_idx)
    }
    # Return final indices
    return(final_idx)
}


extract_data <- function(data, filtered_ids, subsample=500, organism, census, obs_filter, cell_columns, dataset_info, relabel_path,
      dims=20) {

  brain_cell_subsampled_ids = subsample_cells(data, filtered_ids, census, subsample)

          seurat_obj <- get_seurat(
            census = census,
            organism = organism,
            # var_value_filter = gene_filter, # Uncomment if needed
            obs_value_filter = obs_filter,  # Ensure this is constructed correctly
            obs_column_names = cell_columns,
            obs_coords = brain_cell_subsampled_ids,
            var_index="feature_id"
          )
          print("Subsampling successful.")
      # Subset the Seurat object to keep only genes with HGNC"
          # Filter out genes that are not expressed in at least 3 cells
          seurat_obj <- seurat_obj[rowSums(seurat_obj[["RNA"]]@counts > 0) >= 3, ]
          #seurat_obj <- seurat_obj %>% NormalizeData() %>% 
                    #ScaleData() %>% 
                    #FindVariableFeatures() %>% 
                    #RunPCA(npcs=dims) %>% RunUMAP(dims=1:dims)
          #batch_key="dataset_id"
          #seurat_obj <- harmony_integration(seurat_obj, dims=20, integration_key=batch_key)

          newmeta <- left_join(seurat_obj@meta.data , dataset_info, by="dataset_id")  %>%
              select(-soma_joinid.y, everything()) %>%  # Remove soma_joinid.y and keep all other columns
              rename(soma_joinid = soma_joinid.x) 
          rownames(newmeta) <- colnames(seurat_obj)
          seurat_obj@meta.data  <- newmeta
          seurat_obj <- relabel_wrapper(seurat_obj, relabel_path=relabel_path)
          return(seurat_obj)
 }

# Function to extract data based on a specified split value
split_and_extract_data <- function(data, split_column, subsample=500, organism, census, cell_columns, dataset_info, relabel_path, dims=20) {
  # Get unique split values from the specified column
  unique_values <- unique(data[[split_column]])
  # Initialize an empty list to store results
  refs <- list()
  # Loop through each unique split value
  for (split_value in unique_values) {
    # Filter the brain observations based on the split value
    filtered_ids <- data[data[[split_column]] == split_value, ]$soma_joinid
    # Create a filter expression for the Seurat object
    obs_filter <- paste0(split_column, " == '", split_value, "'")
    # Check if the number of filtered ids is greater than or equal to the subsample size
    seurat_obj <- extract_data(data, filtered_ids, subsample, organism, census, obs_filter, cell_columns, dataset_info, relabel_path, dims=dims)
    # Get the unique dataset title
    dataset_titles <- unique(seurat_obj@meta.data[["dataset_title"]])
    # Determine the name to use based on the number of dataset titles
    if (length(dataset_titles) > 1) {
      # Use the split column value as the name
      name_to_use <- split_value
    } else {
      # Use the unique dataset title
      name_to_use <- dataset_titles
    }
    # Store the Seurat object in the refs list using the determined name
    refs[[name_to_use]] <- seurat_obj
  } 
  return(refs)
}

# Call the function with your parameters
get_census <- function(organism="homo_sapiens", census_version="2024-07-01", subsample=500, split_column="dataset_id", dims=20, relabel_path,
                                      ref_collections=c("Transcriptomic cytoarchitecture reveals principles of human neocortex organization",
                                                                                   "SEA-AD: Seattle Alzheimer’s Disease Brain Cell Atlas")) {
      census <- open_soma(census_version)
      # Open obs SOMADataFrame
      dataset_info <- census$get("census_info")$get("datasets")$read()$concat() %>% as.data.frame()
      brain_obs <- census$get("census_data")$get("homo_sapiens")$get("obs")$read(
        value_filter = "tissue_general == 'brain' & is_primary_data == True & disease == 'normal'",
        column_names = c("assay", "cell_type", "sex", "tissue", "tissue_general", "suspension_type", "disease","dataset_id",
        "development_stage","soma_joinid")
      )
      # Concatenates results to an Arrow Table
      brain_obs <-  brain_obs$concat() %>% as.data.frame()
      brain_obs <- left_join(brain_obs, dataset_info, by="dataset_id")
      brain_obs <- brain_obs %>%
        select(-soma_joinid.y, everything()) %>%  # Remove soma_joinid.y and keep all other columns
        rename(soma_joinid = soma_joinid.x) 
      # Allen datasets ?
     # if (organism == "homo_sapiens") {

      brain_obs_filtered <- brain_obs[
          brain_obs$collection_name %in% ref_collections, ]

      brain_obs_filtered <- brain_obs_filtered[!brain_obs_filtered$cell_type %in% c("unknown","glutamatergic neuron"), ]

      #} else if (organism == "mus_musculus") {
      #  brain_obs_filtered <- brain_obs[
      #    brain_obs$collection_name %in% c("x"), ]  # Replace "x" with actual collection names for mouse data
     # } else {
     #   stop("Unsupported organism")
     # }

      # Subsampling cells and extracting the soma_joinid
      set.seed(1)
      # if organism="homo_sapiens set to Homo sapiens"
      #if mus_musculus set to Mus musculus
      if (organism=="homo_sapiens"){
          organism="Homo sapiens"
      } else if (organism=="mus_musculus") {
          organism=="Mus musculus"
      }
      cell_columns <-  c("assay", "cell_type", "tissue", 
      "tissue_general", "suspension_type", 
      "disease","dataset_id","development_stage",
      "soma_joinid")
      refs <- list()
      #everything except whole cortex
      refs <- split_and_extract_data(brain_obs_filtered, split_column=split_column, subsample, organism, census, cell_columns, dataset_info, 
        relabel_path=relabel_path, dims=dims)
      filtered_ids = brain_obs_filtered$soma_joinid
      #whole cortex
      seurat_obj <- extract_data(brain_obs_filtered, filtered_ids, subsample, organism, census, obs_filter=NULL, 
        cell_columns, dataset_info, relabel_path,
        dims=dims) 
      refs[["whole cortex"]] <- seurat_obj

     # names(refs) <- c("whole cortex", unique(brain_obs_filtered$tissue))
    lapply(names(refs), function(x){
      dataset_title <- gsub(" ","_",x)
      #p <- DimPlot(refs[[x]], group.by=c("subclass","assay","tissue","dataset_title"),ncol=2, pt.size=2) +
        #theme(
          #text = element_text(size = 20),  # Adjust the base text size
          #axis.title = element_text(size = 25),  # Adjust axis title size
          #axis.text = element_text(size = 20),  # Adjust axis label size
          #legend.text = element_text(size = 20),  # Adjust legend text size
          #legend.title = element_text(size = 20)  # Adjust legend title size
      #)
      #switch back
      if (organism=="Homo sapiens"){
        organism="homo_sapiens"
      } else if (organism=="Mus musculus") {
        organism=="mus_musculus"
      }
      #ggsave(plot=p, file=file.path(projPath, paste0("/refs/census/",dataset_title,"_",subsample,"_umap.png")), width=35, height=25)
      dir.create("refs", showWarnings = FALSE)
      new_dataset_title <- gsub("[()]", "", gsub("\\/", "_", gsub(" ", "_", gsub("'","",dataset_title))))
      saveRDS(refs[[x]], file=file.path("refs",paste0(new_dataset_title,".rds")))
      #meta <- refs[[x]]@meta.data %>% select("cell_type","class","subclass","family") %>% unique()
      #write.table(meta, file=file.path(projPath, paste0("/meta/relabel/",dataset_title,"_relabel.tsv")), sep="\t",row.names=FALSE)
  })
      
#return(refs)
}