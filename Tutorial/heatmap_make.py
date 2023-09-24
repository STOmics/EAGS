import anndata
import numpy as np
import pandas as pd
import plotly.io as pio
import scanpy as sc
import plotly.graph_objects as go
import os
import scipy.sparse as sp
from scipy.spatial import distance
from sklearn.neighbors import NearestNeighbors
import warnings
from tqdm import tqdm
warnings.filterwarnings("ignore")



def visulize_markers(df:pd.DataFrame,gene:str,save_path:str="./",size:int=2,semititle:str=''):

    color_scale = "inferno"
    if not os.path.isdir(save_path):
        os.makedirs(save_path)
    fig = go.Figure(data = go.Scatter(x = df['x'],
                                  y = df['y'],
                                  mode = 'markers',
                                  opacity=1.0,
                                  marker = dict(
                                  size = size,
                                  color = df[gene],
                                  cmax = 0.9*df[gene].max(),
                                  cmin = 0.1*df[gene].max(),
                                 colorbar=dict(thickness=20, len = 0.3, y = 0.1, tickfont = dict(size = 10), yanchor = "bottom", title = str(semititle)+ str(gene)),
                                 colorscale = color_scale,
                                  showscale = True,
                         reversescale=True
                                  )))
    fig.update_layout(yaxis = dict(title="", linecolor="white", showgrid=False, zeroline = False, visible = False),
                      xaxis = dict(title="", linecolor="white", showgrid=False, zeroline = False, visible = False),
                      autosize=False, width=500, height=380,
                      margin=dict(l=0,r=0,b=0,t=0,pad=0),
                      plot_bgcolor = "rgb(250,250,250)",
                      paper_bgcolor = "rgb(250,250,250)")

    pio.write_image(fig, os.path.join(save_path, str(gene) + "_counts.pdf"))

def show_markerlist(adata:anndata.AnnData,markerlist:list,scale:bool=True,size:int=2,save_path:str="./",semititle:str=''):
    for gene in tqdm(markerlist):
        if gene not in adata.var.index.tolist():
            print(gene)
            continue
        ad = adata[:, gene]
        sc.pp.scale(ad, max_value = 10)
        df = pd.DataFrame(ad.X, columns = [gene], index = ad.obs.index)
        df["x"] = adata.obsm["spatial"][:,0]
        df["y"] = adata.obsm["spatial"][:,1]
        visulize_markers(df=df,gene=gene,save_path=save_path,size=size,semititle=semititle)















