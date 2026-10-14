# Pipeline Defaults vs. Benchmarking: Audit Notes

Updated 2026-09-28. Tables: `BENCHMARKING_AUDIT_TABLES.md`.

## Mouse

- `classifier-effect-glmm`: binomial GLMMs per reference, BH corrected. No outright classifier winner (aggregated: RF 30, Seurat 33, kNN 31 of 94 significant contrasts; Seurat wins outright for only 7/16 subclass types).
- Stable across references: Microglia recall favors Seurat but is unreliable (0 to 100% per study). OPC favors kNN/RF over Seurat (Seurat F1 = 0, aggregated reference only).
- Most other types (Neuron, Glutamatergic, Macrophage, Vascular, Endothelial, Astrocyte, Oligodendrocyte) flip by reference or split between precision and recall, so there is no safe default.
- RF is the most cutoff sensitive classifier, and `cutoff=0` avoids it.
- `subsample_ref=500` is no better than 100, so lower it.
- The cohort is thin (7 studies, 69 samples). Accepted, proceeding.

## Census map fixes

`census_map_mouse_author.tsv` had three silent drop bugs: a raw label with no map row vanishes before subsampling, with no error or warning.

- Microglia (fixed): Micro-PVM was mapped to Macrophage instead of Microglia (isocortex/hippocampal taxonomy only; AUROC .98 to .99 to Microglia in 6/6 query studies). Microglia F1 went from 0.00 to 0.26 before to 0.74 to 0.92 after (whole cortex, subclass, cutoff 0, ref 500).
- OPC (fixed): MOp's `OPC` label (~10k cells) had no row, because the map only had `oligodendrocyte precursor cell` (Tabula Muris spelling). OPC F1 rose for all three classifiers.
- L6 IT Car3 (fixed): MOp's label (720 cells) had no row. It matches isocortex `Car3` at AUROC .9996 to .9999.
- Tabula Muris dropped from `ref_collections`, along with its lowercase ontology name rows. Its only unique labels (Ependymal, Leukocyte) appear in 0/7 query studies.
- Cost: Macrophage F1 falls to near zero for all three classifiers, because its cells are now called Microglia.
- GLMM (q<0.05): Seurat beats kNN on six cell types (recall: Endothelial, GABAergic, Oligodendrocyte, OPC, Vascular; precision: Glutamatergic). kNN beats Seurat on four (Astrocyte precision, Microglia recall, Neural stem cell recall, Neuron). Seurat's edge over RF is larger.
- Sample level, raw (whole cortex, subclass, cutoff 0, ref 100): macro F1 is kNN 0.813, Seurat 0.809, RF 0.740. A plain mean over samples lets studies with more samples count more.
- Sample level, modeled (macro F1 model with study, reference, cutoff, subsample_ref and treatment state terms): kNN 0.887, Seurat 0.860, RF 0.892, with overlapping 95% CIs. RF is lowest in the raw scores but not in the model. The model estimates are marginal means over the other terms, so they are not the same slice as the raw cutoff 0 scores, and the difference is not yet reconciled.
- Seurat is not highest on every metric. Raw (cutoff 0, ref 100) it leads on weighted F1 (0.912 vs 0.901 kNN), micro F1 (0.920 vs 0.912) and accuracy (0.905 vs 0.886), but kNN is slightly ahead on macro F1 (0.813 vs 0.809) and precision (0.954 vs 0.952), and Seurat is lowest in the modeled macro F1 (0.860). Gaps are small and CIs overlap.
- Per cell type recall, precision and F1 by classifier, with significance brackets: `classifier_boxplots/README.md`, and the Seurat vs kNN gaps in F1, recall and precision: `classifier_gap/README.md` (symlink).
- Sweep finished 2026-09-25 on branch `census-map-fixes`. Detail: `cell-type-heterogeneity/FINDINGS.md`.

F1, mean ± sd across queries (whole cortex, subclass, cutoff 0, ref 500):

| Cell type | Classifier | Before | After |
|---|---|---|---|
| Microglia | kNN | 0.00 ± 0.03 | 0.81 ± 0.26 |
| Microglia | RF | 0.00 ± 0.03 | 0.80 ± 0.22 |
| Microglia | Seurat | 0.26 ± 0.36 | 0.92 ± 0.16 |
| OPC | kNN | 0.24 ± 0.31 | 0.79 ± 0.26 |
| OPC | RF | 0.01 ± 0.04 | 0.55 ± 0.35 |
| OPC | Seurat | 0.00 ± 0.03 | 0.82 ± 0.27 |

## Human

RF confirmed best overall default: wins most sig contrasts (35/77) and ties/beats Seurat on F1 for its exception list (precision: L2/3-6 IT, SST, VIP, Glutamatergic-family; recall: deep layer non-IT, L6 CT, L5 ET) - those are small-margin precision/recall losses, not F1 losses. kNN is the clear laggard on F1 (e.g. deep layer non-IT 0.64 vs RF/Seurat ~0.89).

Including Ma et al., decided; actual F1 cost still unmeasured (formula mismatch between compared result sets).

**RESOLVED**: sub500 `whole_cortex` ref cache was stale (zero immune cells); rebuilt 09-15. Verified on clean sub100 data that Seurat's Immune/T-Cell failure is real, not the cache bug: it misclassifies real immune cells as Microglia (88/160) or other neural/glial types, essentially never predicting Immune correctly. scvi_knn/scvi_rf handle it fine. Use `subsample_ref=100` for human GLMM until sub500 is re-scored.

## Open questions

- ~~MOp/10x fixes Microglia (F1 0.07->0.80) but drops OPC + hippocampal taxonomy.~~ RESOLVED 2026-09-22: the OPC drop was the census_map bug above, not a genuine MOp/10x limitation - MOp has ~10k OPC cells under the raw label `OPC`, just never mapped. Fixed; awaiting rerun results.
- Reference breadth vs. accuracy for Gemma uploads: judgment call.
