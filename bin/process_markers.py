#!/user/bin/python3

from pathlib import Path
import os
import sys
import scanpy as sc
import numpy as np
import pandas as pd
import anndata as ad
import re
from scipy.sparse import csr_matrix
import warnings

import adata_functions
from adata_functions import *

import matplotlib.pyplot as plt
import seaborn as sns
import json
import argparse
import os
import json
from types import SimpleNamespace
from scipy.stats import median_abs_deviation
import csv

# Function to parse command line arguments
def parse_arguments():
  parser = argparse.ArgumentParser(description="Classify cells given 1 ref and 1 query")
  parser.add_argument('--markers_file', type=str, default="/space/grp/rschwartz/rschwartz/cell_annotation_cortex.nf/meta/Hybrid_subclass_markers.json")

  if __name__ == "__main__":
      known_args, _ = parser.parse_known_args()
      return known_args

def main():
    
    # Parse command line arguments
  args = parse_arguments()
  markers_file = args.markers_file

  # Load the markers file
  with open(markers_file, 'r') as f:
    markers = json.load(f, object_hook=lambda d: SimpleNamespace(**d))




  with open("/space/grp/rschwartz/rschwartz/cell_annotation_cortex.nf/meta/cell_type_markers.tsv", "w", newline="") as tsvfile:
      writer = csv.writer(tsvfile, delimiter="\t")
      writer.writerow(["Cell Type", "Subtype", "Markers"])  # updated header
      
      for cell_type in markers.cell_types:
          data = cell_type
          if not hasattr(data, 'subtypes'):
              for marker_group in cell_type.markers:
                gene_list = [gene.replace("+", "") for gene in marker_group.genes]
                for gene in gene_list:
                  row = [data.name, subtype.name, gene]
                  writer.writerow(row)

          else:
            for subtype in data.subtypes.cell_types:
              for marker_group in subtype.markers:
                gene_list = [gene.replace("+", "") for gene in marker_group.genes]
                for gene in gene_list:

                  row = [data.name, subtype.name, gene]
                  writer.writerow(row)

      


if __name__ == "__main__":
  main()    # Extract the subclass markers
    