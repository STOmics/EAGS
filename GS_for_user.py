"""
This file stores the original code of the EGAS.
In EGAS you have some choice to your file.
First, you need your file to be h5ad format first.
If not, we have script to exchange gem to h5ad first than using EGAS


This file has been furnished by Lvtongxuan BGI:28650
If you have any Question
please chat with me by using
   Wechat : lvtongxuan666
or Email :  lvtongxuan@genomics.cn

"""
import anndata
import scanpy as sc
import numpy as np
import pandas as pd
from pandas import DataFrame
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.io as pio
import scanpy as sc
import plotly.graph_objects as go
import os
import scipy.sparse as sp
from scipy.spatial import distance
from sklearn.neighbors import NearestNeighbors
import warnings
warnings.filterwarnings("ignore")

def _gaussan_c(dis,gs=0.95,a=1,b=0):
    '''
    :param dis: Smooth weight Distance which has beed calculated
    :param gs: self_adaptively using,that's means the gradient of Gaussian curve (Gaussian parameter c),default value is 0.95
    :param a: default parameter, default value is 1
    :param b: default parameter, default value is 0
    :return: the gradient of Gaussian curve C
    '''
    return np.sqrt(-(dis-b)**2/2/np.log(gs/a))

    
def _gaussan_weight(dis,a=1,b=0,c=12500):
    '''
    :param dis:Smooth weight Distance which has beed calculated
    :param a: default parameter, default value is 1
    :param b: default parameter, default value is 0
    :param c: Control the gradient of Gaussian curve
    :return:
    '''
    gs=a*np.exp(-(dis-b)**2/2/(c**2))
    return gs


def _sp_graph_weight(condition,arr,a=1,b=0,c=12500):
    """
    The diagonal position of the Smooth weight matrix returned by the first
    version of the Smooth matrix may be 0, which may affect subsequent calculations
    :param condition:
    :param arr:
    :param a:default parameter, default value is 1
    :param b:default parameter, default value is 0
    :param c:
    :return:
    """
    out = arr.copy()
    row, col, _=sp.find(condition)
    for ro,co in zip(row, col):
        out[ro,co]=_gaussan_weight(arr[ro,co],a=a,b=b,c=c)
    return out


def _sp_graph_weight_1(condition,arr,dis_filter,a=1,b=0,c=12500):
    """
    Smooth Matrix Second Edition - The diagonal position of the Smooth weight matrix in this version is forcibly set to 1.
    :param condition:
    :param arr:
    :param dis_filter:
    :param a:
    :param b:
    :param c:
    :return:
    """
    out =np.identity(arr.shape[0])
    row, col, _=sp.find(condition)
    for ro,co in zip(row, col):
        if arr[ro,co]>dis_filter:
            out[ro,co]=_gaussan_weight(arr[ro,co],a=a,b=b,c=c)
    return sp.csr_matrix(out)

def gaussian_smooth_fixed(adata:anndata.AnnData,a:float=1,b:float=0,c:float=12500,
                    n_comps:int=50,n_neighbors:int=10,normalize_zscore:bool=True):
    '''
    The curve of the smoothing weight is a fixed value
    :param adata: adata file need to be imputed
    :param a: Gaussian parameter a
    :param b: Gaussian parameter b
    :param c: default value C which control the Gaussian curve
    :param n_comps: dimensions are taken when using principal component analysis
    :param n_neighbors:Calculate the upper limit of neighbors for each cell
    :param normalize_zscore:Whether normalize is required
    :return:Smoothed adata file(h5ad format)
    '''
    adata=adata.copy()

    ### normalize
    if normalize_zscore==True:
        adata.var['mt']=adata.var_names.str.startswith(('MT-','mt-'))
        sc.pp.calculate_qc_metrics(adata,qc_vars=['mt'],inplace=True)
        adata=adata[adata.obs['pct_counts_mt']<10].copy()
        sc.pp.filter_cells(adata,min_counts=300)
        sc.pp.filter_genes(adata,min_cells=10)
        # max_counts=np.percentile(adata.X.sum(1),98.0)
        # sc.pp.filter_cells(adata,max_counts=max_counts)
        # raw_adata=adata.copy().raw.to_adata()
        raw_adata=adata.copy()
        sc.pp.normalize_per_cell(adata,counts_per_cell_after=1e4)
        sc.pp.log1p(adata)
        if type(adata.X)==sp.csr.csr_matrix or type(adata.X)==sp.csc.csc_matrix:
            adata.X = adata.X.toarray()
        adata.X = (adata.X - adata.X.mean(0)) / adata.X.std(0)
    ##### pca
    sc.tl.pca(adata, svd_solver='arpack',n_comps=n_comps)
    Euc_distance=distance.cdist(adata.obsm["spatial"],adata.obsm["spatial"]).astype(np.float32)
    nbrs = NearestNeighbors(n_neighbors=n_neighbors, algorithm='ball_tree').fit(adata.obsm['X_pca'][:,:n_comps])

    adjecent_matrice=nbrs.kneighbors_graph(adata.obsm['X_pca'][:,:n_comps])
    adata.obsp["adjacent_matric"]=adjecent_matrice.astype(np.float32)
    aa=adata.obsp["adjacent_matric"].multiply(Euc_distance) ## 自己和自己的距离为0
    aa=aa.tocsr()

    ##### smoothing
    raw_adata.obsp["gauss_weight"]=_sp_graph_weight(aa>0,aa,a,b,c)
    temp_nor_para=np.squeeze(np.asarray(np.sum(raw_adata.obsp["gauss_weight"],axis=1)))
    new_adata=raw_adata
    new_adata.X=np.asarray(((raw_adata.obsp["gauss_weight"].dot(raw_adata.X)).T/temp_nor_para).T)
    return new_adata


def gaussian_smooth_adaptively(adata:anndata.AnnData,smooth_threshold:float=90,a:float=1,b:float=0,
                    n_comps:int=50,n_neighbors:int=10,normalize_zscore:bool=True):
    """
    The curve of the smoothing weight using adaptive weight calculation
    :param adata:
        adata file need to be imputed
    :param smooth_threshold:
        Control the gradient of Gaussian curve (Gaussian parameter c)
    :param a:
        Gaussian parameter a
    :param b:
        Gaussian parameter b
    :param n_comps:
        Number of principal components to use for calculating neighborhoods
    :param n_neighbors:
        number of nearest neighbors from which to compute kernel bandwidth
    :param normalize_zscore:
        Default preprocessing method for raw counts matrice. If "False", it means input anndata.X data have been normalized.
    :return:
        adata file that has been imputed.
    """
    assert smooth_threshold >=0 and smooth_threshold<=100
    adata=adata.copy()
    
    ### normalize
    if normalize_zscore==True:
        adata.var['mt']=adata.var_names.str.startswith(('MT-','mt-'))
        sc.pp.calculate_qc_metrics(adata,qc_vars=['mt'],inplace=True)
        adata=adata[adata.obs['pct_counts_mt']<10].copy()
        sc.pp.filter_cells(adata,min_counts=300)
        sc.pp.filter_genes(adata,min_cells=10)
        max_counts=np.percentile(adata.X.sum(1),98.0)
        sc.pp.filter_cells(adata,max_counts=max_counts)
        raw_adata=adata.copy().raw.to_adata()
        raw_adata=adata.copy()
        sc.pp.normalize_per_cell(adata,counts_per_cell_after=1e4)
        sc.pp.log1p(adata)
        if type(adata.X)==sp._csr.csr_matrix:
            adata.X = adata.X.toarray()
        adata.X = (adata.X - adata.X.mean(0)) / adata.X.std(0)
    else:
        raw_adata=adata.raw.to_adata()
    
    
    ##### pca
    sc.tl.pca(adata, svd_solver='arpack', n_comps=n_comps)
    Euc_distance=distance.cdist(adata.obsm["spatial"],adata.obsm["spatial"]).astype(np.float32)
    nbrs = NearestNeighbors(n_neighbors=n_neighbors, algorithm='ball_tree').fit(adata.obsm['X_pca'][:,:n_comps])
    adjecent_matrice=nbrs.kneighbors_graph(adata.obsm['X_pca'][:,:n_comps])
    adata.obsp["adjacent_matric"]=adjecent_matrice.astype(np.float32)
    aa=adata.obsp["adjacent_matric"].multiply(Euc_distance) ## 自己和自己的距离为0
    aa=aa.tocsr()  
    aa_nonzero=np.asarray(aa[np.nonzero(aa)])
    dist_threshold=np.percentile(aa_nonzero, smooth_threshold)
    c=_gaussan_c(dist_threshold)
    print('C is :'+str(c))
    ##### smoothing
    raw_adata.obsp["gauss_weight"]=_sp_graph_weight(aa>0,aa,a,b,c)
    temp_nor_para=np.squeeze(np.asarray(np.sum(raw_adata.obsp["gauss_weight"],axis=1)))
    new_adata=raw_adata
    new_adata.X=np.asarray(((raw_adata.obsp["gauss_weight"].dot(raw_adata.X)).T/temp_nor_para).T)
    return new_adata