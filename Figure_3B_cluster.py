#!/usr/bin/env python3

# Converted from Figure_3B_cluster.ipynb.
# Notebook cell order is preserved; Markdown cells are retained as comments.

# %%
import cv2
import stereo as st
from sklearn.cluster import DBSCAN
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context
import os

# %%
data_path = "/mnt/e/stereo/gem/callus_cellbin.selected.gem.gz"

# %%
data = st.io.read_gem(file_path=data_path, bin_type='cell_bins')

# %%
def seperate_samples(data, max_dist=200, SN = None):
    # params: data, stereoExpData contains multiple samples
    # params: max_dist the maximum distance to regard two cells in one sample
    ret = []
    clustering = DBSCAN(eps=200, min_samples=1, algorithm="kd_tree").fit(data.position)
    labels = clustering.labels_
    for x in np.unique(labels):
        cell_list = list(data.cell_names[labels==x])
        data_tmp = data.tl.filter_cells(cell_list=cell_list, inplace=False)
        data_tmp.tl.result['sample'] = [x]*data_tmp.shape[0]
        if type(SN)==str :
            data_tmp.tl.result['SN'] = [SN]*data_tmp.shape[0]
        ret.append(data_tmp)
    return ret

# %%
#这里只有一个芯片，只是以四张芯片举例
data_list = []
SN_list = ["callus_cellbin"]
for i, data_tmp in enumerate([data]):
    data_list.extend(seperate_samples(data_tmp, SN = SN_list[i]))
for data_tmp in data_list:
    print(data_tmp.shape)

# %%
data = st.utils.data_helper.merge(data_list[0],data_list[1],data_list[2],data_list[3],data_list[4],data_list[5],data_list[10],data_list[11],data_list[12],data_list[13],data_list[14],data_list[15],data_list[16],data_list[17])
data.shape

# %%
tmp = np.hstack([data.position, np.array([data.cells.batch]).T])#如果报错，可试试赋予 data.cells.batch 值，类型为列表
tmp = pd.DataFrame(tmp)
tmp.columns = ['x', 'y', 'batch']
tmp.index = data.cell_names

tmp = tmp.astype(np.float64)

# %%
gap_dist = 3500 #图中组织距离
delta_x = 0
for i, x in enumerate(tmp['batch'].unique()):
    mu = tmp.loc[tmp['batch']==x, ['x', 'y']].mean()
    tmp.loc[tmp['batch']==x, ['x', 'y']] -= mu.astype('int')
    tmp.loc[tmp['batch']==x, ['x']] += delta_x
    delta_x+=gap_dist
tmp -= tmp.min()
tmp = tmp.astype(int)
tmp = tmp.loc[data.cell_names, ['x','y']].to_numpy()

# %%
data.position = tmp

# %%
data.tl.cal_qc()
data.plt.violin()
data.plt.spatial_scatter()

# %%
data.tl.raw_checkpoint() ## build a checkpoint to save the original expression data
data.tl.normalize_total()
data.tl.log1p()

# %%
data.tl.pca(use_highly_genes=False, n_pcs=30, res_key='pca')

# %%
data.tl.batches_integrate(pca_res_key='pca', res_key='pca_integrated')

# %%
data.tl.neighbors(pca_res_key='pca_integrated', n_pcs=30, res_key='neighbors_integrated')
data.tl.umap(pca_res_key='pca_integrated', neighbors_res_key='neighbors_integrated', res_key='umap_integrated')
data.plt.batches_umap(res_key='umap_integrated')

# %%
data.tl.spatial_neighbors(neighbors_res_key='neighbors_integrated', res_key='spatial_neighbors')
data.tl.leiden(neighbors_res_key='spatial_neighbors', res_key='spatial_leiden',resolution=2.0)
data.plt.cluster_scatter(res_key='spatial_leiden')

# %%
adata = st.io.stereo_to_anndata(data,flavor='seurat',split_batches=False,output='callus_selected14_pcs30_res2.0_20230925.h5ad')

# %%
adata=sc.read("callus_selected14_pcs30_res2.0_20230925.h5ad")
adata

# %%
adata=sc.read("callus_cellbin_selected14_20230922.h5ad")
adata

# %%
celldf = adata.obs.loc[:,['batch','x','y']]

df_list = []
linedf = celldf[celldf.batch.isin(['0','1','2','3'])]
df_list.append(linedf)

linedf = celldf[celldf.batch.isin(['4','5','6','7'])]
minx = min(linedf.x)
linedf['x'] = linedf['x']-minx
linedf['y'] = linedf['y']+3000
df_list.append(linedf)

linedf = celldf[celldf.batch.isin(['8','9','10','11'])]
minx = min(linedf.x)
linedf['x'] = linedf['x']-minx
linedf['y'] = linedf['y']+6000
df_list.append(linedf)

linedf = celldf[celldf.batch.isin(['12','13'])]
minx = min(linedf.x)
linedf['x'] = linedf['x']-minx
linedf['y'] = linedf['y']+9000
df_list.append(linedf)

celldf = pd.concat(df_list)
celldf = celldf.loc[:,['x','y']]

adata.obs['x'] = celldf['x']
adata.obs['y'] = celldf['y']
adata.obsm['spatial'] = celldf.to_numpy()

# %%
with rc_context({'figure.figsize': (10, 10)}):
    sc.pl.spatial(adata,color='spatial_leiden',spot_size=25,legend_fontsize=15,save="callus_selected14_pcs30_res2.0_spatial_leiden_20230925.pdf")
