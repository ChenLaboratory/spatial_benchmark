#!/usr/bin/env python3
"""
Extract 2D projected cell area (xy footprint) from proseg cell boundary polygons.

proseg-output.zarr stores shapes/cell_boundaries as 2D polygons projected onto
the xy plane.  CRS is cleared before computing area to avoid the geographic-CRS
warning from geopandas.

Output CSV columns: cell, xy_area_um2
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
        print(f"Output already exists, skipping: {args.out}")
        sys.exit(0)

    print(f"Reading zarr: {args.zarr}")
    sdata = spatialdata.read_zarr(args.zarr)

    polys = sdata.shapes["cell_boundaries"]
    polys = polys.set_crs(None, allow_override=True)
    polys["xy_area_um2"] = polys.geometry.area

    out = polys[["cell", "xy_area_um2"]]
    out.to_csv(args.out, index=False)
    print(f"Written: {args.out}  ({len(out)} cells)")


def parse_args():
    parser = argparse.ArgumentParser(description="Compute 2D cell area from proseg boundary polygons.")
    parser.add_argument("-zarr", required=True, help="Path to proseg-output.zarr")
    parser.add_argument("-out",  required=True, help="Output CSV path")
    return parser.parse_args()


if __name__ == "__main__":
    main()
