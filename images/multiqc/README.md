# Single-cell multiQC report

## Outliers
An "outlier" in terms of mitochondrial, ribosomal, or hemoglobin content is >X MADs (median absolute deviation) from the median of pct_counts_mito, pct_counts_ribo, pct_counts_hb. See this single-cell best practices e-book for reference. The number of MADs to threshold at is controlled by the --nmads parameter (default 5).

A "count soutlier" >X MADs from fitted values of a linear model ln(genes + 1) ~ lm(counts + 1)
(see this article for reference). I found this to be helpful with removing likely contamination in dan Goldowitz's data. They use log10, but scanpy default is natural log. I'm not sure if it really matters since data seems to be log-linear.

## Doublets
Doublets are predicted with sc.pp.scrublet.

## Examples
### GSE180670
(add link to benchmarking results)
The reason this dataset performs so poorly in benchmarking is because authors say they sorted for oligodendrocytes. As you can see from the QC report, only about half the cells in every sample are oligodendrocytes. This doesn't seem to be associated with being called an "outlier" on any metric (>5 MADs from median pct_counts_mito|ribo|hb), but that's just an estimate form looking at the report.

I can model the predictive power of these metrics within the benchmarking pipeline if we're interested in that, but of course we can't do that for un-annotated datasets.

I do provide a marker gene heatmap where you can see that a lot of canonical cell type markers for other cells are expressed (it helps if you hit "cluster") on the heatmap. Full disclosure, I got these markers from chatGPT, and while I am familiar with some of them from the literature, I haven't done a deep dive to build a better marker table. It's on my to do list!

### PTSDBrainomics
Example of a "good" dataset according to benchmarking (add link). There are outliers in this dataset as well, I'll need update my benchmarking again to model if they're associated with being classified incorrectly (doubt it)