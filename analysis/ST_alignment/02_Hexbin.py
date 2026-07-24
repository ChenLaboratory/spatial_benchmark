# Bin STalign-aligned iST cells into hexagonal bins matching the Visium spot layout,
# assigning each iST cell a Visium spot barcode.
#
# Usage: python 02_Hexbin.py <sample> <iST_type> <aligned_fname>
#   aligned_fname  output of 01_STalign_iST_to_Visium.py
#                  (<sample>_<iST_type>_STalign_to_Visium_tissue_hires_image_with_point_annotator.csv.gz)

import numpy as np
import pandas as pd
import json
import cv2
from multiprocessing import Pool, cpu_count
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import sys
import itertools
import random
import pickle
from matplotlib.backends.backend_pdf import PdfPages

################ input ####################
sample = sys.argv[1]
iST = sys.argv[2]
aligned_fname = sys.argv[3]


visium_spot_fname = "/spatial/tissue_positions.csv"
visium_json_fname = "/spatial/scalefactors_json.json"

############### output ######################
out_pdf_fname = sample + "_" + iST + "_hexbin_results_overview.pdf"
out_csv_fname = sample + "_" + iST + "_hexbin.csv.gz"

############### functions ####################

def rotate_points(points, alpha, a, b):
    """Rotate points around the point (a, b) by the angle alpha (radians)."""
    rotation_matrix = np.array([
        [np.cos(alpha), -np.sin(alpha)],
        [np.sin(alpha), np.cos(alpha)]
    ])
    translated_points = points - np.array([a, b])
    rotated_points = np.dot(translated_points, rotation_matrix)
    rotated_points += np.array([a, b])
    return np.round(rotated_points, decimals=1)


def hexagon_vertices(center, side_length):
    angle_offset = np.pi / 6  # 30 degrees, for pointy-top hexagons
    angles = np.linspace(0, 2 * np.pi, 6, endpoint=False) + angle_offset
    return np.array([[center[0] + side_length * np.cos(angle), center[1] + side_length * np.sin(angle)] for angle in angles], np.float32)


def point_in_hexagon(point, hexagon):
    return cv2.pointPolygonTest(hexagon, (point[0], point[1]), False) >= 0


def assign_point_to_hexagon(args):
    i, point, hexagon_centers, hexagon_radius = args
    for j, center in hexagon_centers.iterrows():
        hexagon = hexagon_vertices((center["rotated_centre_x"], center["rotated_centre_y"]), hexagon_radius)
        if point_in_hexagon((point["rotated_aligned_x_in_fullres"], point["rotated_aligned_y_in_fullres"]), hexagon):
            return i, center["barcode"]
    return i, ''

#################### determine hexbin centres from Visium spots ####################

visium_spots = pd.read_csv(visium_spot_fname)
visium_spots = visium_spots[visium_spots["in_tissue"] == 1]

visium_x_fullres = np.array(visium_spots['pxl_col_in_fullres'])
visium_y_fullres = np.array(visium_spots['pxl_row_in_fullres'])

with open(visium_json_fname) as json_file:
    json_data = json.load(json_file)
hires_scalef = json_data["tissue_hires_scalef"]

# determine rotation angle from row 0 of the Visium grid
row0 = visium_spots[visium_spots["array_row"] == 0]
row0_y = list(row0["pxl_row_in_fullres"])
row0_x = list(row0["pxl_col_in_fullres"])
phi = np.arctan((row0_y[-1] - row0_y[0]) / (row0_x[-1] - row0_x[0]))

points = np.column_stack((visium_x_fullres, visium_y_fullres))
visium_rotated_points = rotate_points(points, alpha=phi, a=row0_x[0], b=row0_y[0])
visium_rotated_points_int = np.rint(visium_rotated_points).astype(int)
visium_spots["rotated_pxl_col_in_fullres"] = visium_rotated_points_int[:, 0].tolist()
visium_spots["rotated_pxl_row_in_fullres"] = visium_rotated_points_int[:, 1].tolist()

# manually determine spatially equally distributed hexbin centres
col = list(set(visium_spots["array_col"]))
row = list(set(visium_spots["array_row"]))

hexbin_x = np.linspace(min(visium_spots["rotated_pxl_col_in_fullres"]), max(visium_spots["rotated_pxl_col_in_fullres"]), max(col) - min(col) + 1)
hexbin_y = np.linspace(min(visium_spots["rotated_pxl_row_in_fullres"]), max(visium_spots["rotated_pxl_row_in_fullres"]), max(row) - min(row) + 1)

hexbin_x_from_0 = hexbin_x[::2]
hexbin_x_from_1 = hexbin_x[1::2]
hexbin_y_from_0 = hexbin_y[::2]
hexbin_y_from_1 = hexbin_y[1::2]

combinations1 = list(itertools.product(hexbin_x_from_0, hexbin_y_from_0))
combinations2 = list(itertools.product(hexbin_x_from_1, hexbin_y_from_1))
hexbin_centres = pd.DataFrame(combinations1 + combinations2, columns=["rotated_centre_x", "rotated_centre_y"])

# hexbin radius/side length, from spacing of Visium spot grid
v_dis = hexbin_y[1] - hexbin_y[0]
h_dis = (hexbin_x[1] - hexbin_x[0]) * 2
v_dis_R = v_dis * 2 / 3
h_dis_R = h_dis / np.sqrt(3)
hexagon_side_length = min(v_dis_R, h_dis_R)

# subset manual centres having a corresponding Visium spot inside
subset_centres = []
for i, center in hexbin_centres.iterrows():
    hexagon = hexagon_vertices((center["rotated_centre_x"], center["rotated_centre_y"]), hexagon_side_length)
    for j, spot in visium_spots.iterrows():
        if point_in_hexagon((spot["rotated_pxl_col_in_fullres"], spot["rotated_pxl_row_in_fullres"]), hexagon):
            subset_centres.append({
                "index": j,
                "rotated_centre_x": center["rotated_centre_x"],
                "rotated_centre_y": center["rotated_centre_y"],
                "barcode": spot["barcode"]
            })
            break
subset_centres = pd.DataFrame(subset_centres)

with open('subset_centres.pkl', 'wb') as f:
    pickle.dump(subset_centres, f)

############## process iST aligned data ##############

iST_cells = pd.read_csv(aligned_fname)

# convert aligned coordinates (hires-image pixels) into full-resolution pixels
iST_cells["aligned_x_in_fullres"] = iST_cells["aligned_x"] / hires_scalef
iST_cells["aligned_y_in_fullres"] = iST_cells["aligned_y"] / hires_scalef

# rotate aligned coordinates of iST with the same angle as Visium
points = iST_cells[["aligned_x_in_fullres", "aligned_y_in_fullres"]]
iST_rotated_points = rotate_points(points, alpha=phi, a=row0_x[0], b=row0_y[0])
iST_cells["rotated_aligned_x_in_fullres"] = iST_rotated_points[:, 0].tolist()
iST_cells["rotated_aligned_y_in_fullres"] = iST_rotated_points[:, 1].tolist()

############## assign iST cells to hexbins (Visium spots) ##############

args = [(i, point, subset_centres, hexagon_side_length) for i, point in iST_cells.iterrows()]
iST_cells["Visium_spot_id"] = [''] * len(iST_cells)

num_cores = min(40, cpu_count())
with Pool(processes=num_cores) as pool:
    results = pool.map(assign_point_to_hexagon, args)

for i, hex_id in results:
    iST_cells.at[i, "Visium_spot_id"] = hex_id

###################### export results #######################

iST_cells.to_csv(out_csv_fname, compression='gzip')

with PdfPages(out_pdf_fname) as pdf:

    fig, ax = plt.subplots(figsize=(20, 20))
    ax.scatter(np.array(subset_centres["rotated_centre_x"]), np.array(subset_centres["rotated_centre_y"]),
               c='black', s=2, label="manually set hexbin centres")
    ax.scatter(np.array(visium_spots["rotated_pxl_col_in_fullres"]), np.array(visium_spots["rotated_pxl_row_in_fullres"]),
               c='red', s=2, label="Visium: rotated spots under fullres")
    lgnd = plt.legend(loc="upper right", scatterpoints=1, fontsize=10)
    for handle in lgnd.legend_handles:
        handle.set_sizes([10.0])
    ax.set_aspect('equal')
    plt.title('Manually calculated hexbin centres', fontsize=40)
    pdf.savefig(fig)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(20, 20))
    for j, center in subset_centres.iterrows():
        hexagon = hexagon_vertices((center["rotated_centre_x"], center["rotated_centre_y"]), hexagon_side_length)
        ax.add_patch(patches.Polygon(hexagon, closed=True, edgecolor='r', facecolor='none'))
    ax.scatter(np.array(visium_spots["rotated_pxl_col_in_fullres"]), np.array(visium_spots["rotated_pxl_row_in_fullres"]),
               c='black', s=2, label="Visium: rotated spots under fullres")
    lgnd = plt.legend(loc="upper right", scatterpoints=1, fontsize=10)
    for handle in lgnd.legend_handles:
        handle.set_sizes([10.0])
    ax.set_aspect('equal')
    plt.title('Hexbin pattern', fontsize=40)
    pdf.savefig(fig)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(20, 20), dpi=100)
    ax.scatter(np.array(visium_spots["rotated_pxl_col_in_fullres"]), np.array(visium_spots["rotated_pxl_row_in_fullres"]),
               c='red', label="Visium: rotated spots under fullres")
    ax.scatter(iST_rotated_points[:, 0].tolist(), iST_rotated_points[:, 1].tolist(),
               s=1, alpha=0.2, color='blue', label="iST: rotated alignment under fullres")
    lgnd = plt.legend(loc="upper right", scatterpoints=1, fontsize=10)
    for handle in lgnd.legend_handles:
        handle.set_sizes([10.0])
    plt.title('iST aligned to Visium (fullres, rotated)', fontsize=40)
    pdf.savefig(fig)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(40, 40))
    for j, center in subset_centres.iterrows():
        hexagon = hexagon_vertices((center["rotated_centre_x"], center["rotated_centre_y"]), hexagon_side_length)
        ax.add_patch(patches.Polygon(hexagon, closed=True, edgecolor='r', facecolor='none'))
    ax.scatter(np.array(visium_spots["rotated_pxl_col_in_fullres"]), np.array(visium_spots["rotated_pxl_row_in_fullres"]),
               c='black', s=2, label="Visium: rotated spots under fullres")

    iST_cells_plot = iST_cells[iST_cells["Visium_spot_id"] == ""]
    for i, point in iST_cells_plot.iterrows():
        ax.plot(point["rotated_aligned_x_in_fullres"], point["rotated_aligned_y_in_fullres"], 'o', color="blue", markersize=1)
    ax.plot([], [], 'o', color="blue", markersize=5, label="iST cells without hexbin id")

    iST_cells_plot = iST_cells[iST_cells["Visium_spot_id"].isin(random.sample(list(visium_spots["barcode"]), 30))]
    for i, point in iST_cells_plot.iterrows():
        ax.plot(point["rotated_aligned_x_in_fullres"], point["rotated_aligned_y_in_fullres"], 'o', color="red", markersize=1)
    ax.plot([], [], 'o', color="red", markersize=5, label="iST cells inside randomly selected hexbins")

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, loc='upper right', scatterpoints=1, fontsize=10)
    ax.set_aspect('equal')
    plt.title('Hexagon Binning of Points', fontsize=40)
    pdf.savefig(fig)
    plt.close(fig)
