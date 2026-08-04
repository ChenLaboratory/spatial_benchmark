# Align Xenium/MERSCOPE (iST) cell centroids onto the matched Visium H&E image using STalign.
# Landmark points (<sample>_<iST_type>_points.npy, <sample>_Visium_points.npy) are pre-selected
# manually (e.g. with STalign's point_annotator.py) and are a required input to this script.

import numpy as np
import pandas as pd
import torch
import pickle
from STalign import STalign
import matplotlib.pyplot as plt

iST_type = "Xenium"  # or "MERSCOPE"
sample_id = "sample_id"

# iST --------------------------------------------------------------------

if iST_type == "Xenium":
    iST_fname = "/cells.csv.gz"
    df1 = pd.read_csv(iST_fname)
    xI = np.array(df1['x_centroid'])
    yI = np.array(df1['y_centroid'])
elif iST_type == "MERSCOPE":
    iST_fname = "/cell_metadata.csv.gz"
    df1 = pd.read_csv(iST_fname)
    xI = np.array(df1['center_x'])
    yI = np.array(df1['center_y'])
else:
    raise ValueError("invalid iST type")

XI, YI, I, fig = STalign.rasterize(xI, yI, dx=2)
I = np.vstack((I, I, I))  # make into 3xNxM
I = STalign.normalize(I)

np.savez(sample_id + "_" + iST_type, x=XI, y=YI, I=I)

# Visium -------------------------------------------------------------------

Visium_image_file = '/spatial/tissue_hires_image.png'
V = plt.imread(Visium_image_file)
Vnorm = STalign.normalize(V)
J = Vnorm.transpose(2, 0, 1)

YJ = np.array(range(J.shape[1])) * 1.
XJ = np.array(range(J.shape[2])) * 1.

np.savez(sample_id + '_Visium', x=XJ, y=YJ, I=J)

# landmark points ------------------------------------------------------------

pointsIlist = np.load(sample_id + "_" + iST_type + '_points.npy', allow_pickle=True).tolist()
pointsJlist = np.load(sample_id + '_Visium_points.npy', allow_pickle=True).tolist()

pointsI = []
pointsJ = []
for i in pointsIlist.keys():
    for j in range(len(pointsIlist[i])):
        pointsI.append([pointsIlist[i][j][1], pointsIlist[i][j][0]])
for i in pointsJlist.keys():
    for j in range(len(pointsJlist[i])):
        pointsJ.append([pointsJlist[i][j][1], pointsJlist[i][j][0]])

pointsI = np.array(pointsI)
pointsJ = np.array(pointsJ)

# landmark-based LDDMM alignment ---------------------------------------------

if torch.cuda.is_available():
    torch.set_default_device('cuda:0')
    device = 'cuda:0'
else:
    torch.set_default_device('cpu')
    device = 'cpu'

# compute initial affine transformation from points
L, T = STalign.L_T_from_points(pointsI, pointsJ)
A = STalign.to_A(torch.tensor(L), torch.tensor(T))

# keep all other parameters default
params = {'L': L, 'T': T,
          'niter': 200,
          'pointsI': pointsI,
          'pointsJ': pointsJ,
          'device': device,
          'sigmaP': 2e-1,
          'sigmaM': 0.18,
          'sigmaB': 0.18,
          'sigmaA': 0.18,
          'diffeo_start': 100,
          'epL': 5e-11,
          'epT': 5e-4,
          'epV': 5e1
          }

out = STalign.LDDMM([YI, XI], I, [YJ, XJ], J, **params)

with open('out_objs_using_landmark_2um.pkl', 'wb') as f:
    pickle.dump(out, f)

A = out['A']
v = out['v']
xv = out['xv']

# apply transform to original cell coordinates
tpointsI = STalign.transform_points_source_to_target(xv, v, A, np.stack([yI, xI], 1))
if tpointsI.is_cuda:
    tpointsI = tpointsI.cpu()

# switch from row column coordinates (y,x) to (x,y)
xI_LDDMM = tpointsI[:, 1]
yI_LDDMM = tpointsI[:, 0]

# export aligned coordinates --------------------------------------------------

df3 = pd.DataFrame({
    "aligned_x": xI_LDDMM.cpu() if tpointsI.is_cuda else xI_LDDMM,
    "aligned_y": yI_LDDMM.cpu() if tpointsI.is_cuda else yI_LDDMM,
})
results = pd.concat([df1, df3], axis=1)

results.to_csv(sample_id + "_" + iST_type + "_STalign_to_Visium_tissue_hires_image_with_point_annotator.csv.gz",
               compression='gzip')
