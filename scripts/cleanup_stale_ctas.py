#!/usr/bin/env python3
"""Manual maintenance tool: delete stale-version pipeline cell type assignments
(CTAs) for a Gemma study. Not invoked by the pipeline itself - run by hand when
a repeatedly-reannotated study (e.g. a dev/test fixture) has accumulated CTAs
from old pipeline versions that were never replaced, since gemma-cli's
-replaceCta only replaces an exact name match and this pipeline's CTA names
embed the version (sc-pipeline-<version>-<level>).

Usage:
    python scripts/cleanup_stale_ctas.py --study_name DevBrain --version 2.0.0
    python scripts/cleanup_stale_ctas.py --study_name DevBrain --version 2.0.0 --dry_run
"""
import argparse
import subprocess
import gemmapy


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--study_name", required=True, help="Gemma study short name or ID")
    p.add_argument("--version", required=True, help="Current pipeline version to keep (e.g. 2.0.0)")
    p.add_argument("--protocol_prefix", default="sc-pipeline")
    p.add_argument("--use_staging", action="store_true", help="Use Gemma staging instead of production")
    p.add_argument("--dry_run", action="store_true", help="List what would be deleted without deleting")
    return p.parse_args()


def main():
    args = parse_args()
    gp = gemmapy.GemmaPy(path="staging" if args.use_staging else "production")

    try:
        dim = gp.raw.get_dataset_single_cell_dimension(dataset=args.study_name).to_dict()["data"]
    except Exception as e:
        print(f"No existing single-cell dimension for {args.study_name}, nothing to clean up: {e}")
        return

    gemma_cmd = "gemma-cli-staging" if args.use_staging else "gemma-cli"

    # Match on our prefix, not a specific level: anything named "<prefix>-<x>" that
    # isn't in the CURRENT version's namespace is stale, whatever level "<x>" ends
    # in - so a level this script has never heard of still gets cleaned up.
    current_prefix = f"{args.protocol_prefix}-{args.version}-"
    own_prefix = f"{args.protocol_prefix}-"

    for cta in dim["cell_type_assignments"]:
        name = cta["name"]
        if not name.startswith(own_prefix) or name.startswith(current_prefix):
            continue
        if args.dry_run:
            print(f"[dry run] would delete: {name} (id={cta['id']})")
            continue
        print(f"Deleting stale CTA: {name} (id={cta['id']})")
        subprocess.run(
            [gemma_cmd, "deleteSingleCellData", "-e", args.study_name, "-deleteCta", str(cta["id"])],
            check=True,
        )


if __name__ == "__main__":
    main()
