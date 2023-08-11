"""
This test function has been furnished by Lvtongxuan BGI:28650
If you have any Question
please chat with me by using
   Wechat : lvtongxuan666
or Email :  lvtongxuan@genomics.cn
This fuction help you to use Gaussian-smoothing without choose parameter and normalize adata.h5ad by yourself

if you want to make it better
please check the README.md and set parameter by yourself
"""

import os
import sys
import scanpy as sc
import numpy as np
import scanpy
import anndata
import scipy.sparse as sp
import warnings
warnings.filterwarnings("ignore")
import scanpy as sc
import numpy as np
import pandas as pd
from GS_for_user import gaussian_smooth_fixed,gaussian_smooth_adaptively
import scipy.sparse
import argparse

if __name__ == '__main__':
    ## add input adata path(h5ad format)
    adata_path = "your adata file path"
    ## adata_file path
    adata = sc.read_h5ad(adata_path)
    ## fixed smooth weight EGAS
    adata_after_smooth = gaussian_smooth_adaptively(adata=adata,
                                                    smooth_threshold=90,
                                                    a=1,
                                                    b=0,
                                                    n_comps=50,
                                                    n_neighbors=10,
                                                    normalize_zscore=True)
    ## save the result
    adata_after_smooth.write_h5ad('output_file_path')
