# Harmonic analysis of spatial transcriptomics data
[![Biorxiv badge](https://zenodo.org/badge/doi/10.1101/2023.06.30.547258.svg)](https://doi.org/10.1101/2023.06.30.547258) ⬅️ manuscript <br>

This repository includes
1. a minimal package for performing filtering
2. all [notebooks](notebooks) used to generate figures

See instructions below for installing and using the package.

## 1. Installation
```
pip install git+https://github.com/wanglab-broad/harmonics@main
```

## 2. Usage
This package can be used in Python as follows:
```python
import scanpy as sc
import harmonics as hm

# Load data
adata = sc.read_h5ad('/path/to/your/data.h5ad')

# High-pass filter genes using a square-root kernel
hm.sfilter(
    adata,
    'juxtacrine',
)
```
For a list of parameters, their meanings, and their types, run `help(hm.sfilter)`.

Alternatively, one can run `hm.sfilter` from the shell:
```bash
harmonics \
--read_path /path/to/your/data.h5ad \
--write_path /path/to/your/filtered_data.h5ad \
--filter_keys juxtacrine \
```
For a list of parameters in the shell, run `harmonics --help`.
