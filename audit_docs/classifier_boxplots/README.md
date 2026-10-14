# Classifier comparison by cell type (subclass)

Mouse census map fixes sweep, whole cortex reference, cutoff 0, `subsample_ref` 100. Regenerated from `celltype-differences` with the plotting scripts sorted by support.

![Subclass recall](boxplots_subclass_recall.png)

![Subclass precision](boxplots_subclass_precision.png)

![Subclass F1](f1_paired_subclass.png)

**What they are.** Per cell type recall, precision and F1 for scVI kNN, scVI RF and Seurat. Recall and precision come from `celltype-differences/glmm-precision-recall` (`boxplots_subclass_*.png`), F1 from `celltype-differences/paired-f1-wilcoxon` (`f1_paired_subclass.png`).

**How to read them.** Rows are cell types grouped by lineage and sorted by total support within each lineage, largest at the top (cells summed over queries, counted once per query, shown in each label). The three classifiers are dodged inside each row. Large points are means, bars are mean ± SD (recall and precision figures), and small points are per sample values. The F1 figure also lists the number of queries in each label. Brackets compare each classifier pair: `*` means q < 0.05 and `ns` means not significant. The recall and precision brackets come from a binomial GLMM with a study random effect (cells weighted). The F1 brackets come from a paired Wilcoxon on per query F1 (each query weighted equally, study ignored). Missing brackets mean the pair could not be tested.

**What they show.** Recall and F1 agree in direction for most cell types, but not always in significance. Microglia kNN vs Seurat is significant for recall (GLMM) and not for F1 (Wilcoxon), because the large studies dominate the cell weighted GLMM while two small studies favor Seurat.
