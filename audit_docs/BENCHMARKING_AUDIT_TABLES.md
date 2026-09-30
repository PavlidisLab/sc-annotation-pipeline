# Benchmarking audit: supporting tables

Appendix to `BENCHMARKING_AUDIT.md`, which is prose-only. By-study detail
that `classifier-effect-glmm` doesn't output (it works on pooled counts, not
per-study breakdowns), backing specific claims in the Microglia and Neural
stem cell sections there. The raw per-reference F1 tables previously here
were removed once `classifier-effect-glmm` covered all 4 mouse references
with real significance testing, superseding a plain-F1 comparison.

## Ground-truth cell counts (Microglia, OPC, Neural stem cell)

Same across all methods and references, since these are true label counts,
not predictions:

| study | Microglia | OPC | Neural stem cell |
|---|---|---|---|
| GSE124952 | 76 | 32 | 0 |
| GSE181021.2 | 21 | 57 | 0 |
| GSE185454 | 21 | 3 | 0 |
| GSE214244.1 | 21 | 12 | 0 |
| GSE247339.1 | 206 | 232 | 15 |
| GSE247339.2 | 289 | 220 | 81 |
| **total** | **634** | **556** | **96** |

## Microglia recall by study and reference

Seurat's Microglia recall, by study and reference (10x and SMART-Seq v4
carry no Microglia class at all, so Macrophage dominates there by default,
an absence, not a failure):

| study | TBI? | MOp atlas | `aggregated` |
|---|---|---|---|
| GSE124952 | no | 97% | 36% |
| GSE181021.2 | no | 100% | 95% |
| GSE185454 | no | 0% | 81% |
| GSE214244.1 | no | 0% | 43% |
| GSE247339.1 | yes | 71% | 7% |
| GSE247339.2 | yes | 67% | 4% |

## Neural stem cell precision/recall/F1 by reference and method

Unlike Microglia and OPC, where specific classifiers just never predict the
label, all three classifiers get real, nonzero recall for Neural stem cell
wherever it's a class at all:

| reference | method | GSE247339.1 (precision / recall / F1) | GSE247339.2 (precision / recall / F1) |
|---|---|---|---|
| MOp atlas | scvi_knn | 14% / 80% / 0.23 | 50% / 84% / 0.63 |
| MOp atlas | scvi_rf | 7% / 67% / 0.13 | 41% / 90% / 0.56 |
| MOp atlas | seurat | 4% / 33% / 0.06 | 23% / 38% / 0.29 |
| `aggregated` | scvi_knn | 19% / 40% / 0.26 | 64% / 68% / 0.66 |
| `aggregated` | scvi_rf | 4% / 53% / 0.07 | 17% / 33% / 0.23 |
| `aggregated` | seurat | 17% / 33% / 0.23 | 57% / 47% / 0.51 |
