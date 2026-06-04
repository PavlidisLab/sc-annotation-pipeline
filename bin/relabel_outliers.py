#!/usr/bin/env python
"""Relabel QC-outlier cells as "unknown" in a cell type assignment (CTA) file.

Joins a per-study outlier mask onto a classifier CTA file and overwrites the
predicted label (and ontology URI) with "unknown" for any cell flagged as an
outlier in any QC category. This reuses the same "unknown" label the classifier
already assigns to below-cutoff cells, so downstream consumers see a single
"not confidently called" category whether the reason was low probability or QC.

Safe to apply post-hoc: the annotation path is per-cell independent (the scVI
encoder is frozen and the random forest predicts per cell), so relabeling an
outlier never changes any other cell's label.
"""

import argparse

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--celltype_file', type=str, required=True,
                        help='CTA TSV with columns sample_id, cell_id, cell_type, cell_type_uri')
    parser.add_argument('--mask_file', type=str, required=True,
                        help='Mask TSV with columns sample_id, cell_id, category, value (true/false)')
    parser.add_argument('--output', type=str, required=True,
                        help='Output path for the relabeled CTA TSV')
    return parser.parse_args()


def main():
    args = parse_args()

    celltypes = pd.read_csv(args.celltype_file, sep='\t')
    mask = pd.read_csv(args.mask_file, sep='\t')

    # Identify outlier cells (value == "true" in any category; the mask file already
    # collapses all QC categories into a single "mask" row per cell).
    outliers = mask.loc[mask['value'].astype(str).str.lower() == 'true', ['sample_id', 'cell_id']].copy()

    # Join on sample_id + cell_id using string keys so an int/float dtype mismatch
    # (e.g. 483963 vs 483963.0) can't silently drop the join.
    for col in ['sample_id', 'cell_id']:
        celltypes['_' + col] = celltypes[col].astype(str)
        outliers['_' + col] = outliers[col].astype(str)
    outlier_keys = set(zip(outliers['_sample_id'], outliers['_cell_id']))

    is_outlier = celltypes.apply(
        lambda r: (r['_sample_id'], r['_cell_id']) in outlier_keys, axis=1
    )

    celltypes.loc[is_outlier, 'cell_type'] = 'unknown'
    if 'cell_type_uri' in celltypes.columns:
        celltypes.loc[is_outlier, 'cell_type_uri'] = pd.NA

    celltypes = celltypes.drop(columns=['_sample_id', '_cell_id'])
    celltypes.to_csv(args.output, sep='\t', index=False)

    print(f"Relabeled {int(is_outlier.sum())} / {len(celltypes)} cells as 'unknown' (QC outliers).")


if __name__ == '__main__':
    main()
