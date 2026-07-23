#!/usr/bin/env python3

import argparse
from census_utils import setup


def parse_arguments():
    parser = argparse.ArgumentParser(description="Download scVI model for a given organism and census version.")
    parser.add_argument('--organism', type=str, default='homo_sapiens')
    parser.add_argument('--census_version', type=str, default='2024-07-01')
    return parser.parse_args()


args = parse_arguments()
model_path = setup(organism=args.organism, version=args.census_version)
