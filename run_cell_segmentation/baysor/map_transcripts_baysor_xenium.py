#!/usr/bin/env python3
"""
Convert Baysor segmentation.csv (Xenium) to 10X MTX format for Seurat/Scanpy.

Xenium-specific behaviour:
  - Unassigned filter : cell == 0  (integer sentinel; NaN/empty not used)
  - Extra quality filter : qv >= 20  (Phred-scaled detection quality; Xenium-specific column)
  - Cell IDs in barcodes.tsv : prefixed with 'cell_'  (e.g. cell_12345)
  - Control feature types :
      NegControlProbe_ / antisense_  -> "Negative Control Probe"
      NegControlCodeword_            -> "Negative Control Codeword"
      BLANK_                         -> "Blank Codeword"
"""

import argparse
import csv
import os
import sys
import numpy as np
import pandas as pd
import scipy.sparse as sparse
import scipy.io as sio
import subprocess


def main():
    # Parse input arguments.
    args = parse_args()

    # Check for existence of input file.
    if not os.path.exists(args.baysor):
        print(f"The specified Baysor output ({args.baysor}) does not exist!")
        sys.exit(1)

    # Check if output folder already exists.
    if os.path.exists(args.out):
        print(f"The specified output folder ({args.out}) already exists!")
        sys.exit(1)

    # Read required columns
    print("Reading transcripts...")
    transcripts_df = pd.read_csv(args.baysor,
                                 usecols=["gene", "cell", "assignment_confidence", "qv"])

    print(f"Filtering transcripts (qv >= 20, confidence >= {args.conf_cutoff})...")

    # Xenium: unassigned transcripts have cell == 0 (integer sentinel).
    # qv (Phred-scaled quality value) >= 20 is the standard Xenium threshold for reliable detections.
    transcripts_df = transcripts_df[
        transcripts_df["cell"].notna() &
        (transcripts_df["cell"] != 0) &
        (transcripts_df["cell"] != "") &
        (transcripts_df["qv"] >= 20) &
        (transcripts_df["assignment_confidence"] >= args.conf_cutoff)
    ].copy()

    print("Aggregating counts...")
    # Group by gene and cell to get the counts (replaces the manual loop)
    # This creates a long-format table of [gene, cell, count]
    counts = transcripts_df.groupby(['gene', 'cell']).size().reset_index(name='count')

    # Convert to categorical to map names to integer indices for the sparse matrix
    counts['gene'] = pd.Categorical(counts['gene'])
    counts['cell'] = pd.Categorical(counts['cell'])

    features = counts['gene'].cat.categories
    cells = counts['cell'].cat.categories

    # Build the sparse matrix directly using coordinate (COO) format
    # (row_index, col_index, value)
    row_ind = counts['gene'].cat.codes
    col_ind = counts['cell'].cat.codes
    data = counts['count']

    sparse_mat = sparse.coo_matrix((data, (row_ind, col_ind)),
                                    shape=(len(features), len(cells)))

    # Write output files
    write_sparse_mtx(args, sparse_mat, cells, features)
    print(f"Successfully wrote matrix to {args.out}")


#--------------------------
# Helper functions

def parse_args():
    """Parses command-line options for main()."""
    summary = 'Map Xenium/Baysor transcripts to a Seurat/Scanpy-compatible matrix.'
    parser = argparse.ArgumentParser(description=summary)

    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument('-baysor', required=True,
                               help="Path to the *segmentation.csv file.")
    requiredNamed.add_argument('-out', required=True,
                               help="Name of output folder.")

    parser.add_argument('-conf_cutoff', default=0, type=float,
                        help="Confidence threshold (default: 0)")

    # rep_int is kept for backwards compatibility, though no longer needed for the new logic
    parser.add_argument('-rep_int', default=100000, type=int, help="Reporting interval.")

    return parser.parse_args()


def write_sparse_mtx(args, sparse_mat, cells, features):
    """Write matrix in Seurat/Scanpy-compatible MTX format"""

    if not os.path.exists(args.out):
        os.mkdir(args.out)

    # Write matrix in MTX format.
    sio.mmwrite(os.path.join(args.out, "matrix.mtx"), sparse_mat)

    # Write cells as barcodes.tsv.
    with open(os.path.join(args.out, "barcodes.tsv"), 'w', newline='') as tsvfile:
        writer = csv.writer(tsvfile, delimiter='\t', lineterminator='\n')
        for cell in cells:
            # Match the prefixing you used in your original script
            writer.writerow(["cell_" + str(cell)])

    # Write features as features.tsv.
    with open(os.path.join(args.out, "features.tsv"), 'w', newline='') as tsvfile:
        writer = csv.writer(tsvfile, delimiter='\t', lineterminator='\n')
        for f in features:
            feature = str(f)
            if feature.startswith("NegControlProbe_") or feature.startswith("antisense_"):
                writer.writerow([feature, feature, "Negative Control Probe"])
            elif feature.startswith("NegControlCodeword_"):
                writer.writerow([feature, feature, "Negative Control Codeword"])
            elif feature.startswith("BLANK_"):
                writer.writerow([feature, feature, "Blank Codeword"])
            else:
                writer.writerow([feature, feature, "Gene Expression"])

    # Gzip files for Seurat compatibility
    print("Compressing files...")
    subprocess.run(f"gzip -f {args.out}/*", shell=True)


if __name__ == "__main__":
    main()
