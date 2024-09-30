# TODO: support for user-defined filter kernels

"""
Perform filtering of gene expression signals over a tissue graph domain.
"""

from typing import Optional, Collection

import scanpy as sc
import numpy as np
from scipy.spatial import Delaunay
import scipy.sparse as sp
import pygsp as pg
from tqdm import tqdm


def clean(
    adata: sc.AnnData,
    min_genes: int,
    min_counts: int,
    min_cells: int,
    percent_cells: float,
):
    """\
    Clean the data by filtering cell/genes out based on a minimum gene/cell thresholds.

    Parameters
    ----------
    adata
        AnnData object
    min_genes
        The minimum number of genes each cell must express
    min_counts
        The minimum number of counts each cell must express
    min_cells
        The minimum number of cells each gene must be expressed in
    percent_cells
        Sets `min_cells` based on a percentage of the total cells in the dataset
    
    Returns
    -------
    Modifies `adata` by removing thresholded cells/genes in place.
    """

    if min_genes:
        print(f'\tRemoving cells with fewer than {min_genes} genes', flush=True)
        sc.pp.filter_cells(adata, min_genes=min_genes)
    if min_counts:
        print(f'\tRemoving cells with fewer than {min_counts} counts', flush=True)
        sc.pp.filter_cells(adata, min_counts=min_counts)
    if min_cells:
        print(f'\tRemoving genes found in fewer than {min_cells} cells', flush=True)
        sc.pp.filter_genes(adata, min_cells=min_cells)
    if percent_cells:
        print(f'\tRemoving genes found in less than {percent_cells} percent of cells', flush=True)
        sc.pp.filter_genes(adata, min_cells=len(adata)*percent_cells)


def normalize(
    adata: sc.AnnData,
    normalize_total: bool,
    log1p: bool,
    scale: bool,
    layer: str,
):
    """\
    Normalize the data.

    Parameters
    ----------
    adata
        AnnData object
    normalize_total
        Whether to normalize the counts per cell
    log1p
        Whether to log transform the normalized counts per gene
    scale
        Whether to scale the log transformed expression per gene
    layer 
        Which layer in `adata.layers` to use (`None` corresponds to adata.X)
    
    Returns
    -------
    Modifies `adata` by normalizing expression values in place.
    """

    if normalize_total:
        print('\tnormalize_total', flush=True)
        sc.pp.normalize_total(adata, layer=layer)
    if log1p:
        print('\tlog1p', flush=True)
        sc.pp.log1p(adata, layer=layer)
    if scale:
        print('\tscale', flush=True)
        sc.pp.scale(adata, layer=layer)


def delaunay(
    coords: np.ndarray,
) -> sp.csr_matrix:
    """\
    Compute a Delaunay triangulation over a set of coordinates (i.e. tissue mesh).

    Parameters
    ----------
    coords
        An n x d numpy array corresponding to n cells in d spatial dimensions
    
    Returns
    -------
    A
        Sparse adjacency matrix corresponding to the triangulation
    """

    # Adapted from Squidpy's squidpy.gr._build under BSD-3-Clause License
    # Source: https://github.com/scverse/squidpy/blob/main/src/squidpy/gr/_build.py
    tri = Delaunay(coords)
    n_cells = len(coords)
    indptr, indices = tri.vertex_neighbor_vertices
    A = sp.csr_matrix(
        (np.ones_like(indices, dtype=np.float64), indices, indptr),
        shape=(n_cells, n_cells),
    )

    return A


def _iter_filter(
    X: np.ndarray,
    g: pg.filters.Filter,
    sparse: bool,
    order: int,
    pbar: bool,
) -> np.ndarray:
    """\
    Iteratively filter gene expression signals over the graph, one gene at a time.
    Saves memory (at the cost of time) compared to filtering the entire matrix at once.

    Parameters
    ----------
    X 
        The gene expression matrix
    g
        A PyGSP filter of type `pygsp.filters.Filter`
    sparse
        Whether the gene expression matrix `X` is sparse
    order
        The order of the approximation used by PyGSP (i.e. of the Chebyshev polynomial)
    pbar
        Whether to display a progress bar corresponding to the number of genes filtered
    
    Returns
    -------
    X_filtered
        A filtered version of `X`.
        Necessarily dense, as the result of filtering will inherently not be sparse.
    """

    X_filtered = np.zeros(X.shape)
    n_genes = X.shape[1]
    if pbar:
        iterator = tqdm(range(n_genes))
    else:
        iterator = range(n_genes)
    for i in iterator:
        x = X[:,i]
        if sparse:
            x = np.array(x.todense()).flatten()
        X_filtered[:,i] = g.filter(x, order=order)

    return X_filtered


def sfilter(
    adata: sc.AnnData,
    filter_keys: Collection[str],
    layer: Optional[str] = None,
    sparse: Optional[bool] = False,
    lap_type: Optional[str] = 'normalized',
    order: Optional[int] = 30,
    spatial_key: Optional[str] = 'spatial',
    tau: Optional[float] = 5,
    pbar: Optional[bool] = True,
):
    """\
    Filter gene expression signals over a spatial domain.

    Parameters
    ----------
    adata 
        AnnData object
    filter_keys
        Types of filtering desired
        Options are
            - 'heat': low-pass filter used for region identification
            - 'juxtacrine': high-pass filter used for interaction identification
            - 'juxtacrine_dual': low-pass dual filter of 'juxtacrine'
            - 'paracrine': mid-pass filter used for boundary identification
    layer
        Layer in `adata.layers` to filter
    sparse
        Whether the matrix of signals (e.g. `adata.X`) is sparse
    lap_type
        Laplacian type to use (i.e. 'combinatorial' or 'normalized')
    order
        The order of the approximation used by PyGSP (i.e. of the Chebyshev polynomial)
    spatial_key
        The field under `adata.obsm` under which the spatial coordinates are stored
    tau
        Parameter determining the extent of diffusion
    pbar
        Whether to display a progress bar corresponding to the number of genes filtered
    
    Returns
    -------
    Modifies adata in place by adding filtered expression signals under
    `adata.layers[filter_key]` for each key in `filter_keys`.
    """

    # Additional input parsing
    if isinstance(filter_keys, str):
        filter_keys = [filter_keys]

    # Choose signals to filter
    if layer:
        X = adata.layers[layer]
    else:
        X = adata.X

    # Create domain
    print('Creating graph domain', flush=True)
    G = pg.graphs.Graph(delaunay(adata.obsm[spatial_key]), lap_type=lap_type)
    G.estimate_lmax()

    # Create filter(s)
    gdict = dict()
    if 'heat' in filter_keys:
        gdict['heat'] = pg.filters.Filter(G, lambda x: np.exp(-tau*x))
    if 'juxtacrine' in filter_keys:
        gdict['juxtacrine'] = pg.filters.Filter(G, lambda x: np.sqrt(x))
    if 'juxtacrine_dual' in filter_keys:
        gdict['juxtacrine_dual'] = pg.filters.Filter(G, lambda x: np.sqrt(G.lmax-x))
    if 'paracrine' in filter_keys:
        gdict['paracrine'] = pg.filters.Filter(G, lambda x: np.exp(-tau*x) * np.sqrt(x))

    # Apply filter(s) to signals
    for filter_key in filter_keys:
        g = gdict[filter_key]
        print(f'Filtering using a {filter_key} kernel', flush=True)
        adata.layers[filter_key] = _iter_filter(X, g, sparse, order, pbar)
