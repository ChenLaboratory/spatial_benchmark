import matplotlib.pyplot as plt
import os
from scipy.stats import mode
import sys
import pandas as pd
import numpy as np
import matplotlib.colors as mcolors
import cv2
import argparse

def load_cmap(cmap_tsv):
    """
    Load cmap dictionary from TSV:
    cluster_label   color
  
    """
    df = pd.read_csv(cmap_tsv, sep="\t")
    color_map_dict={}
    for n,x,y,z in zip(df['Name'],df['R'],df['G'],df['B']):
        color_map_dict[str(n)] = (mcolors.to_hex((x,y,z)))
    print(color_map_dict)
    return color_map_dict




def plot_pixel_cli(_args):

    parser = argparse.ArgumentParser(prog="plot_pixel")

    # Input / Output
    parser.add_argument("--input_tsv", required=True,
                        help="Input TSV containing X,Y,cluster labels")
    parser.add_argument("--output", required=True,
                        help="Output image path (PNG)")

    # Columns
    parser.add_argument("--x_col", default="X",
                        help="Column name for X coordinates")
    parser.add_argument("--y_col", default="Y",
                        help="Column name for Y coordinates")
    parser.add_argument("--cluster_col", default="K1",
                        help="Column name for cluster labels")

    # Plot settings
    parser.add_argument("--pixel_size", type=float, default=10,
                        help="Grid size for pixelation")
    parser.add_argument("--min_points", type=int, default=5,
                        help="Minimum points per pixel to color")

    parser.add_argument("--join_method", choices=["avg", "major"],
                        default="major",
                        help="How to color pixels")

    parser.add_argument("--background_color", default="black",
                        help="Background color")

    # Colormap + Highlighting
    parser.add_argument("--cmap_dict_tsv", required=True,
                        help="TSV mapping cluster → hex color")

    parser.add_argument("--highlight_clusters", nargs="*",
                        help="Clusters to highlight (others become grey)")

    args = parser.parse_args(_args)

    # -------------------------
    # Load input + cmap
    # -------------------------
    tx_meta = pd.read_csv(args.input_tsv, sep="\t",skiprows=3)  
    cmap = load_cmap(args.cmap_dict_tsv)
    print(cmap)
    tx_meta[args.cluster_col] = tx_meta[args.cluster_col].astype(str)  # Ensure cluster labels are strings for consistent mapping

    # Convert highlight list
    highlight_clusters = None
    if args.highlight_clusters:
        highlight_clusters = list(set(args.highlight_clusters))
    print('highlight_clusters',highlight_clusters)
    # -------------------------
    # Compute pixel indices
    # -------------------------
    tx_meta["px"] = (tx_meta[args.x_col] // args.pixel_size).astype(int)
    tx_meta["py"] = (tx_meta[args.y_col] // args.pixel_size).astype(int)

    # -------------------------
    # Pixel coloring function
    # -------------------------
    def pixel_color(df):

        if len(df) < args.min_points:
            return pd.Series([np.nan, np.nan, np.nan])

        labels = df[args.cluster_col]

        # Highlight mode
        if highlight_clusters is not None:

            # Split points into highlighted vs non-highlighted
            is_highlight = labels.astype(str).isin(highlight_clusters)

            # If no highlight points → pale grey fallback
            if is_highlight.sum() == 0:
                return pd.Series(mcolors.to_rgb("#7d7d7d91"))

            # Assign weights
            weights = np.where(is_highlight, 1.0, 0.05)

            # Convert cluster labels → RGB
            rgb_vals = labels.map(
                lambda lbl: mcolors.to_rgb(
                    cmap.get(lbl, args.background_color)
                )
            )

            rgb_array = np.vstack(rgb_vals)

            # Weighted mean
            weighted_rgb = np.average(rgb_array, axis=0, weights=weights)

            return pd.Series(weighted_rgb)

        # Default behavior (no highlighting)
        if args.join_method == "major":
            major_label = labels.mode().iloc[0]
            return pd.Series(
                mcolors.to_rgb(cmap.get(major_label, args.background_color))
            )

        else:  # avg
            rgb_vals = labels.map(
                lambda lbl: mcolors.to_rgb(
                    cmap.get(lbl, args.background_color)
                )
            )
            rgb_array = np.vstack(rgb_vals)
            return pd.Series(rgb_array.mean(axis=0))
    # -------------------------
    # Group by pixel + build image
    # -------------------------
    pixel_df = (
        tx_meta.groupby(["px", "py"])
        .apply(pixel_color)
        .reset_index()
    )

    pixel_df.columns = ["px", "py", "R", "G", "B"]

    # Pivot into image grid
    pivot_r = pixel_df.pivot(index="py", columns="px", values="R").iloc[::-1]
    pivot_g = pixel_df.pivot(index="py", columns="px", values="G").iloc[::-1]
    pivot_b = pixel_df.pivot(index="py", columns="px", values="B").iloc[::-1]

    img = np.dstack([pivot_r, pivot_g, pivot_b])

    # Fill NaNs with background
    bg_rgb = np.array(mcolors.to_rgb(args.background_color))
    nan_mask = np.isnan(img)
    img[nan_mask] = np.take(bg_rgb, nan_mask.nonzero()[2])

    # Save image
    img_uint8 = np.clip(img * 255, 0, 255).astype(np.uint8)
    img_bgr = cv2.cvtColor(img_uint8, cv2.COLOR_RGB2BGR)

    cv2.imwrite(args.output, img_bgr)
    print(f"✅ Saved pixelated image → {args.output}")


if __name__ == "__main__":
    plot_pixel_cli(sys.argv[1:])
