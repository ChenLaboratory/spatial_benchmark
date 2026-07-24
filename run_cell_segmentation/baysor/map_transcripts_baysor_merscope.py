#!/usr/bin/env python3
"""
Convert Baysor segmentation.csv (MERSCOPE) to 10X MTX format for Seurat/Scanpy.

MERSCOPE-specific behaviour:
  - Unassigned filter : cell is NaN or empty string  (no integer sentinel like Xenium)
  - No extra quality filter  (no qv column in MERSCOPE Baysor output)
  - Cell IDs in barcodes.tsv : written as-is  (alphanumeric, e.g. CR6d88aca73-1)
  - Control feature types :
      Blank-*  -> "Blank Codeword"
      (not specifying NegControlProbe / NegControlCodeword categories in MERSCOPE)
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
    args = parse_args()

    if not os.path.exists(args.baysor):
        print(f"The specified Baysor output ({args.baysor}) does not exist!")
        sys.exit(1)

    if os.path.exists(args.out):
        print(f"The specified output folder ({args.out}) already exists!")
        sys.exit(1)

    # Only load the three columns needed; avoids reading large unused fields (x, y, fov, etc.)
    print("Reading transcripts...")
    transcripts_df = pd.read_csv(args.baysor,
                                 usecols=["gene", "cell", "assignment_confidence"])

    print(f"Filtering transcripts (confidence >= {args.conf_cutoff})...")

    # MERSCOPE: unassigned transcripts have NaN/empty cell (Xenium uses integer sentinel 0).
    # Low-confidence assignments are excluded to reduce noise in the count matrix.
    transcripts_df = transcripts_df[
        transcripts_df["cell"].notna() &
        (transcripts_df["cell"] != "") &
        (transcripts_df["assignment_confidence"] >= args.conf_cutoff)
    ].copy()

    print("Aggregating counts...")
    # Count transcripts per (gene, cell) pair — equivalent to a UMI count matrix.
    counts = transcripts_df.groupby(["gene", "cell"]).size().reset_index(name="count")

    # Convert gene/cell names to integer category codes so they can index the sparse matrix.
    counts["gene"] = pd.Categorical(counts["gene"])
    counts["cell"] = pd.Categorical(counts["cell"])

    features = counts["gene"].cat.categories  # ordered gene names → row labels
    cells    = counts["cell"].cat.categories  # ordered cell IDs  → column labels

    # Build sparse matrix in COO format: (row=gene index, col=cell index, value=count).
    # COO is efficient to construct from triples; scipy converts to CSC for mmwrite.
    row_ind = counts["gene"].cat.codes
    col_ind = counts["cell"].cat.codes
    data    = counts["count"]

    sparse_mat = sparse.coo_matrix((data, (row_ind, col_ind)),
                                    shape=(len(features), len(cells)))

    write_sparse_mtx(args, sparse_mat, cells, features)
    print(f"Successfully wrote matrix to {args.out}")


def parse_args():
    summary = "Map MERSCOPE/Baysor transcripts to a Seurat/Scanpy-compatible matrix."
    parser = argparse.ArgumentParser(description=summary)

    req = parser.add_argument_group("required named arguments")
    req.add_argument("-baysor", required=True,
                     help="Path to the Baysor segmentation.csv file.")
    req.add_argument("-out", required=True,
                     help="Name of output folder.")

    parser.add_argument("-conf_cutoff", default=0, type=float,
                        help="Confidence threshold (default: 0)")

    return parser.parse_args()


def write_sparse_mtx(args, sparse_mat, cells, features):
    os.mkdir(args.out)

    # matrix.mtx: sparse count matrix in Matrix Market format (genes × cells)
    sio.mmwrite(os.path.join(args.out, "matrix.mtx"), sparse_mat)

    # barcodes.tsv: one cell ID per line.
    # MERSCOPE cell IDs are alphanumeric strings (e.g. CR6d88aca73-1); written as-is.
    # Xenium uses integer IDs that need a 'cell_' prefix — not needed here.
    with open(os.path.join(args.out, "barcodes.tsv"), "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        for cell in cells:
            writer.writerow([str(cell)])

    # features.tsv: gene_id, gene_name, feature_type (3-column format expected by Read10X).
    # MERSCOPE blank codewords are named 'Blank-*'; all other entries are gene expression.
    # Xenium has additional types (NegControlProbe_, NegControlCodeword_, BLANK_) not present here.
    with open(os.path.join(args.out, "features.tsv"), "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        for feat in features:
            feat_str = str(feat)
            if feat_str.startswith("Blank-"):
                writer.writerow([feat_str, feat_str, "Blank Codeword"])
            else:
                writer.writerow([feat_str, feat_str, "Gene Expression"])

    # Gzip all three files — required for Seurat's Read10X() to recognise them
    print("Compressing files...")
    subprocess.run(f"gzip -f {args.out}/*", shell=True)


if __name__ == "__main__":
    main()
