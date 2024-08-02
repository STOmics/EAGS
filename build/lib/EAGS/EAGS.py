"""
This test function has been furnished by Lvtongxuan BGI:28650
If you have any Question
please chat with me by using
   Wechat : lvtongxuan666
or Email :  lvtongxuan@genomics.cn

if you want to make it better
please check the README.md and set parameter by yourself
"""

import warnings

from EAGS.EAGS_function import gaussian_smooth_adaptively

warnings.filterwarnings("ignore")
import scanpy as sc
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='EAGS using')
    parser.add_argument('--input', type=str, required=True,
                        default='/home/users/input.h5ad',
                        help="input Anndata file name")
    parser.add_argument('--smooth_threshold', type=float, default=90, help="Control the gradient of Gaussian curve")
    parser.add_argument('--a', type=float, default=1, help='parameter: a(default 1) in Gaussian weight')
    parser.add_argument('--b', type=float, default=0, help='parameter: b(default 0) in Gaussian weight')
    parser.add_argument('--n_comps', type=int, default=50, help='number of principal components to use for calculating neighborhoods')
    parser.add_argument('--n_neighbors', type=int, default=10, help='number of nearest neighbors from which to compute kernel bandwidth')
    parser.add_argument('--normalize_zscore', type=bool, default=False, help='default preprocessing method for raw counts matrice')
    parser.add_argument('--output', type=str, default='/home/users/gaussian.h5ad', help='output Anndata file name (defaut ./gauss.ha5ad)')
    args = parser.parse_args()
    adata = sc.read_h5ad(args.input)
    gaussian = gaussian_smooth_adaptively(adata=adata,
                                          smooth_threshold=args.smooth_threshold,
                                          a=args.a,
                                          b=args.b,
                                          n_comps=args.n_comps,
                                          n_neighbors=args.n_neighbors,
                                          normalize_zscore=args.normalize_zscore)
    gaussian.write_h5ad(args.output)
