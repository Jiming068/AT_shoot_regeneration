#!/usr/bin/env python3

# Converted from Figure_S7B_in_situ_expression.ipynb.
# Notebook cell order is preserved; Markdown cells are retained as comments.

# %%
import scanpy as sc
import pandas as pd
import numpy as np
#import squidpy as sq
import scanpy.external as sce
import matplotlib.pyplot as plt
import scvelo as scv
import os

# %%
adata = sc.read('/data/work/callus_cellbin_selected14_subtype_ms.h5ad')

# %%
tissue_figsize = {'D6-root-3-1':(5,5),
'D6-root-1-3':(4.6,4.2),
'D8-root-2-1':(3.9,3.2),
'D6-root-4-0':(6.5,4.5),
'D6-root-1-2':(4.5,4.8),
'D6-root-4-3':(5.8,4.8),
'D10-root-1-1':(4.7,5.4),
'D8-root-2-0':(4.2,2.6),
'D6-root-3-2':(4.9,5.3),
'D8-root-1-1':(4.1,3.5),
'D10-root-1-0':(5,5.2),
'D6-root-1-1':(4.9,4),
'D6-root-1-4':(4.3,5.1),
'D8-root-1-0':(4.3,3.5)}

# %%
import matplotlib as mpl
from matplotlib.colors import ListedColormap, LinearSegmentedColormap, to_rgb, cnames, is_color_like
cmaptest = LinearSegmentedColormap.from_list("mycmap", ['#fff7f3ff','#8e0152ff'])
cmaptest

# %%
set(adata.obs['tissue'])

# %%
set_type = '20260402_somatic_embryo'

# %%
genedf1 = pd.read_csv('/data/work/somatic_embryo_genes.txt', sep ='\t')

# %%
for tissue in set(adata.obs['tissue']):
    if not os.path.exists('./figures/'+set_type+'/'+tissue+'/'):
        os.makedirs('./figures/'+set_type+'/'+tissue+'/', exist_ok=True)
    try:
        msadata = sc.read('fromSAP/'+tissue+'.Ms.h5ad')
    except OSError as e:
        print(f"Skip {tissue} due to file error: {e}")
        continue  # 跳过该组织，继续下一个

    msadata.obs['celltype'] = adata[msadata.obs_names].obs['celltype']
    for num in list(genedf1.index):
        gene_id = genedf1.loc[num,'ID']
        file_name = './figures/'+set_type+'/'+tissue+'/'+str(num)+'_'+'_'.join(list(genedf1.loc[num,['Symbol','ID']]))
        try:
            scv.pl.scatter(msadata, basis='spatial', color=[gene_id], size=50, dpi=600, figsize=tissue_figsize[tissue], alpha=0.9,
                color_map=['#fff7f3ff','#8e0152ff'], layer='Ms', smooth=True, save=file_name+'.Ms_smooth09.pdf', show=False, colorbar=False)
        except:
            print(gene_id)

# %%
import zipfile
import os

# %%
def zip_folder(folder_path, output_zip):
    with zipfile.ZipFile(output_zip, 'w', zipfile.ZIP_DEFLATED) as zipf:
        for root, dirs, files in os.walk(folder_path):
            for file in files:
                file_path = os.path.join(root, file)
                # 在 zip 中保存相对路径，使压缩包解压后保持目录结构
                arcname = os.path.relpath(file_path, start=folder_path)
                zipf.write(file_path, arcname)

# %%
# 使用
zip_folder('/data/work/figures', '/data/work/figures.zip')
