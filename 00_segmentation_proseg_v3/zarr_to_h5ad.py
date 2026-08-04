#!/usr/bin/env python3
"""
Convert proseg zarr output to h5ad for downstream Seurat conversion.
"""

import argparse
import os
import sys
import spatialdata


def main():
    args = parse_args()

    if not os.path.exists(args.zarr):
        print(f"Zarr store not found: {args.zarr}")
        sys.exit(1)

    if os.path.exists(args.out):
        print(f"Output file already exists: {args.out}")
        sys.exit(1)

    print(f"Reading zarr: {args.zarr}")
    sdata = spatialdata.read_zarr(args.zarr)

    print(f"Writing h5ad: {args.out}")
    sdata.tables["table"].write_h5ad(args.out)
    print("Done.")


def parse_args():
    parser = argparse.ArgumentParser(description="Convert proseg zarr to h5ad.")
    parser.add_argument("-zarr", required=True, help="Path to proseg-output.zarr")
    parser.add_argument("-out",  required=True, help="Output h5ad file path")
    return parser.parse_args()


if __name__ == "__main__":
    main()
