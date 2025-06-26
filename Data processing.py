import random
import numpy as np
import anndata as ad
import scanpy as sc
import scanpy.external as sce
import muon as mu
from muon import atac as ac
from matplotlib import pyplot as plt

random.seed(1)
np.random.seed(1)

rna_raw = ad.read("data/0920_mergeRNA.h5ad")
rna_raw

atac_raw = ad.read("data/0920_mergeATAC.h5ad")
atac_raw

mdata = mu.MuData({'rna': rna_raw, 'atac': atac_raw})
print(mdata.var_names[mdata.var_names.str.startswith('RPS') + mdata.var_names.str.startswith('RPL')])
mdata = mdata[:, ~(mdata.var_names.str.startswith('RPS') + mdata.var_names.str.startswith('RPL'))]
mdata.var_names_make_unique()
mdata

mdata.update()
mu.pp.intersect_obs(mdata)
mdata

## RNA
rna = mdata.mod['rna']

# Normalisation
rna.layers["counts"] = rna.X.copy()
sc.pp.normalize_total(rna, target_sum=1e4)
sc.pp.log1p(rna)
rna.layers["lognorm"] = rna.X.copy()

# Define informative features
sc.pp.highly_variable_genes(rna, min_mean=0.02, max_mean=4, min_disp=0.5)
print(f"the number of highly variable genes: {np.sum(rna.var.highly_variable)}")

# Scaling
rna.raw = rna
sc.pp.scale(rna, max_value=10)

## ATAC
atac = mdata.mod['atac']

# Normalisation
atac.layers["counts"] = atac.X.copy()
sc.pp.normalize_total(atac, target_sum=1e4)
sc.pp.log1p(atac)
atac.layers["lognorm"] = atac.X.copy()

# Define informative features
sc.pp.highly_variable_genes(atac, min_mean=0.05, max_mean=1.5, min_disp=.5)
print(f"the number of highly variable peaks: {np.sum(atac.var.highly_variable)}")

# Scaling
atac.raw = atac
sc.pp.scale(atac, max_value=10)

# Multi-omics integration
mdata.update()
mu.pp.intersect_obs(mdata)
mdata

# mofa
mu.tl.mofa(mdata, outfile="results/merge_mofa_model_filtered.hdf5", gpu_mode=False, seed=1)

sce.pp.harmony_integrate(mdata, key="sample", basis='X_mofa')
sc.pp.neighbors(mdata, use_rep="X_pca_harmony", random_state=1)
sc.tl.umap(mdata, random_state=1)
mdata.obsm["X_mofa_umap"] = mdata.obsm["X_umap"]
sc.pl.umap(mdata, color="sample", size=2)

sc.tl.leiden(mdata, key_added='leiden_joint', resolution=.7, random_state=1)
sc.pl.umap(mdata, color="leiden_joint", legend_loc='on data', size=2)

mdata.write("results/merge_afterMOFA_filtered.h5mu")



