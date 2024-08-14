"""
Command-line interface.
"""

from argparse import ArgumentParser
from typing import Collection

import scanpy as sc

from .filtering import clean, normalize, sfilter


def main(
    read_path: str,
    write_path: str,
    filter_keys: Collection[str],
    layer: str,
    min_genes: int,
    min_counts: int,
    min_cells: int,
    percent_cells: float,
    sparse: bool,
    normalize_total: bool,
    log1p: bool,
    scale: bool,
    lap_type: str,
    order: int,
    spatial_key: str,
    tau: float,
    pbar: bool,
):
    """\
    Main function for submitting filtering jobs from the shell.
    """

    # Read data
    print(f'Reading {read_path}', flush=True)
    adata = sc.read_h5ad(read_path)

    # Clean data
    if min_genes or min_counts or min_cells or percent_cells:
        print('Cleaning', flush=True)
        print('\tBefore: ', f'{adata.shape[0]} cells x {adata.shape[1]} genes')
        clean(adata, min_genes, min_counts, min_cells, percent_cells)
        print('\tAfter: ', f'{adata.shape[0]} cells x {adata.shape[1]} genes')

    # Normalize data
    if normalize_total or log1p or scale:
        print('Normalizing')
        normalize(adata, normalize_total, log1p, scale, layer)

    # Spatially filter data
    sfilter(
        adata=adata,
        filter_keys=filter_keys,
        layer=layer,
        sparse=sparse,
        lap_type=lap_type,
        order=order,
        spatial_key=spatial_key,
        tau=tau,
        pbar=pbar,
    )

    # Write data
    adata.write(write_path)
    print(f'Written to {write_path}', flush=True)


parser = ArgumentParser()
parser.add_argument('--read_path', type=str)
parser.add_argument('--write_path', type=str)
parser.add_argument('--filter_keys', type=str, nargs='+')
parser.add_argument('--layer', type=str, default=None)
parser.add_argument('--min_genes', type=int, default=None)
parser.add_argument('--min_counts', type=int, default=1)
parser.add_argument('--min_cells', type=int, default=None)
parser.add_argument('--percent_cells', type=float, default=None)
parser.add_argument('--sparse', action='store_true')
parser.add_argument('--normalize_total', action='store_true')
parser.add_argument('--log1p', action='store_true')
parser.add_argument('--scale', action='store_true')
parser.add_argument('--lap_type', type=str, default='combinatorial')
parser.add_argument('--order', type=int, default=30)
parser.add_argument('--spatial_key', type=str, default='spatial')
parser.add_argument('--tau', type=float, default=5)
parser.add_argument('--pbar', action='store_true')
args = parser.parse_args()

main(
    args.read_path,
    args.write_path,
    args.filter_keys,
    args.layer,
    args.min_genes,
    args.min_counts,
    args.min_cells,
    args.percent_cells,
    args.sparse,
    args.normalize_total,
    args.log1p,
    args.scale,
    args.lap_type,
    args.order,
    args.spatial_key,
    args.tau,
    args.pbar,
)
