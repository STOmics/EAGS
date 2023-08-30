# EAGS: efficient and adaptive Gaussian smoothing applied to high-resolved spatial transcriptomics
## Introduction
EAGS is a smoothing approach for spatial transcriptome data with ultrahigh resolution. 
Through the principal component information of the profiling of gene expression and the 
personally identifiable information between cells, it uses a two-factor adaptive smoothing 
weight mechanism to determine the neighborhood relationship of cells. Then it uses the 
self-developed smoothing weight formula to calculate the smoothing contribution and 
recalculate the gene expression of each cell.

This work belongs to the thematic series of Spatio-temporal omics algorithms.
## Schematic Diagram
![img.png](img/fig1.jpg)

Figure: A. The spatio-temporal omics data processing flow. B. The EAGS method flow.


## Quick Start

### Download the GitHub Repository
[Download](https://github.com/STOmics/EAGS/archive/refs/heads/main.zip) this GitHub repository, and extract the contents into a folder.


### Install Dependencies
```bash
### Python enviroment constructed by Conda
conda create -n EAGS python=3.8.13
conda activate EAGS
pip install -r requirements.txt
```

### Running EAGS Script from the Command-line

```bash
cd ./EAGS-main/

# EAGS
python EAGS.py
--input input.h5ad
--smooth_threshold 90
--a 1
--b 0
--n_comps 50
--n_neighbors 10
--normalize_zscore True
--output gaussian.h5ad
```

EAGS can also be performed directly without some parameter setting as:
```bash
# EAGS
cd ./EAGS-main/

python EAGS.py
--input input.h5ad
--output gaussian.h5ad
```

### Parameter Explanation

The core function of the EAGS method is the `gaussian_smooth_adaptively` function provided in `EAGS.py`, and the 
detailed formula of this function has been given in the `EAGS_function.py`.

The following code is a quick example of running our Smoothing strategy. The function `gaussian_smooth_adaptively()` takes 
in `anndata(H5AD format)` 
object with the cell information(such as gene MID Counts in `anndata.X` , spatial coordinates in `anndata.obsm['Spatial']`). 
For more details on the input anndata(H5AD format), please check on [link](https://anndata.readthedocs.io/en/latest/).

    adata_after_smooth = gaussian_smooth_adaptively(adata=adata,
                                                    smooth_threshold=90,
                                                    a=1,
                                                    b=0,
                                                    n_comps=50,
                                                    n_neighbors=10,
                                                    normalize_zscore=True)

The parameters of  `gaussian_smooth_adaptively` are:
- `adata`: Adata file (H5AD format) need to be smoothed.
- `smooth_threshold`: Control the gradient of Gaussian curve.
- `a`: The a-value of the Gaussian kernel function, which used to control the height of the Gaussian curve.
- `b`: The b-value of the Gaussian kernel function, which used to control the centre of the Gaussian curve.
- `n_comps`: Number of principal components to use for calculating neighborhoods. 
- `n_neighbors`: Number of nearest neighbors from which to compute kernel bandwidth.
- `normalize_zscore`: Default preprocessing method for raw counts matrice. If "False", it means input anndata.X data 
have been normalized.

The output of `gaussian_smooth_adaptively()` is a `anndata(H5AD format)` includes:
- `adata.X`: MID counts have been smoothed.
- `adata.obsm[spatial]`: Spot's spatial coordinates.


## Contact
Any questions or suggestions on EAGS are welcomed! Please report it on issues, 
or contact Tongxuan Lv (lvtongxuan@genomics.cn) or Ying Zhang (zhangying7@genomics.cn).
We recommend using [STOmics Cloud Platform](https://cloud.stomics.tech/) to access and use it.
