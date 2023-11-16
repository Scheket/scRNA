#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 10 12:29:41 2023

@author: liza
"""
import numpy as np
from scipy.sparse import csc_matrix
import scanpy as sc
import pandas as pd

#%% 1st type of files from seurat object.

# THE PARAMETERS of the object is the following:
    # nFeature_RNA > 200 & 
                  # nFeature_RNA < 4000 &
                  # nCount_RNA > 200 &
                  # nCount_RNA < 23000 &
                  # percent.mt < 13 &
                  # percent.ribo < 35

# Read the .mtx file along with the barcodes and gene names
adata = sc.read_mtx("/home/liza/Documents/PhD/scRNAseq_Inbal/All.mtx")
genes = pd.read_csv("/home/liza/Documents/PhD/scRNAseq_Inbal/genes.tsv", header=None, sep='\t')[0]
barcodes = pd.read_csv("/home/liza/Documents/PhD/scRNAseq_Inbal/barcodes.tsv", header=None, sep='\t')[0]

adata = adata.transpose()

# Check dimensions and assign names
if adata.shape[1] == len(genes):
    # If the number of columns in adata.X matches the number of genes, assign var_names
    adata.var_names = genes
    adata.obs_names = barcodes
else:
    # If there's a mismatch, raise an error or alert
    raise ValueError("Mismatch between the number of genes and the number of columns in the AnnData object.")


# Ensure that var_names and obs_names do not have duplicates
adata.var_names_make_unique(join="-")
adata.obs_names_make_unique(join="-")


#%% FROM SEURAT FILES GC

# Read the .mtx file along with the barcodes and gene names
adata_GC = sc.read_mtx("/home/liza/Documents/PhD/scRNAseq_Inbal/GC.mtx")
genes = pd.read_csv("/home/liza/Documents/PhD/scRNAseq_Inbal/GCgenes.tsv", header=None, sep='\t')[0]
barcodes = pd.read_csv("/home/liza/Documents/PhD/scRNAseq_Inbal/GCbarcodes.tsv", header=None, sep='\t')[0]

adata_GC = adata_GC.transpose()

# Check dimensions and assign names
if adata_GC.shape[1] == len(genes):
    # If the number of columns in adata.X matches the number of genes, assign var_names
    adata_GC.var_names = genes
    adata_GC.obs_names = barcodes
else:
    # If there's a mismatch, raise an error or alert
    raise ValueError("Mismatch between the number of genes and the number of columns in the AnnData object.")


# Assign the gene names and cell barcodes
adata_GC.var_names = genes
adata_GC.obs_names = barcodes

# Ensure that var_names and obs_names do not have duplicates
adata_GC.var_names_make_unique(join="-")
adata_GC.obs_names_make_unique(join="-")

#%% FROM SEURAT FILES GP

# Read the .mtx file along with the barcodes and gene names
adata_GP = sc.read_mtx("/home/liza/Documents/PhD/scRNAseq_Inbal/GP.mtx")
genes = pd.read_csv("/home/liza/Documents/PhD/scRNAseq_Inbal/GPgenes.tsv", header=None, sep='\t')[0]
barcodes = pd.read_csv("/home/liza/Documents/PhD/scRNAseq_Inbal/GPbarcodes.tsv", header=None, sep='\t')[0]

adata_GP = adata_GP.transpose()

# Check dimensions and assign names
if adata_GP.shape[1] == len(genes):
    # If the number of columns in adata.X matches the number of genes, assign var_names
    adata_GP.var_names = genes
    adata_GP.obs_names = barcodes
else:
    # If there's a mismatch, raise an error or alert
    raise ValueError("Mismatch between the number of genes and the number of columns in the AnnData object.")


# Assign the gene names and cell barcodes
adata_GP.var_names = genes
adata_GP.obs_names = barcodes

# Ensure that var_names and obs_names do not have duplicates
adata_GP.var_names_make_unique(join="-")
adata_GP.obs_names_make_unique(join="-")
#%%FROM SEURAT FILES WP

# Read the .mtx file along with the barcodes and gene names
adata_WP = sc.read_mtx("/home/liza/Documents/PhD/scRNAseq_Inbal/WP.mtx")
genes = pd.read_csv("/home/liza/Documents/PhD/scRNAseq_Inbal/WPgenes.tsv", header=None, sep='\t')[0]
barcodes = pd.read_csv("/home/liza/Documents/PhD/scRNAseq_Inbal/WPbarcodes.tsv", header=None, sep='\t')[0]

adata_WP = adata_WP.transpose()

# Check dimensions and assign names
if adata_WP.shape[1] == len(genes):
    # If the number of columns in adata.X matches the number of genes, assign var_names
    adata_WP.var_names = genes
    adata_WP.obs_names = barcodes
else:
    # If there's a mismatch, raise an error or alert
    raise ValueError("Mismatch between the number of genes and the number of columns in the AnnData object.")


# Assign the gene names and cell barcodes
adata_WP.var_names = genes
adata_WP.obs_names = barcodes

# Ensure that var_names and obs_names do not have duplicates
adata_WP.var_names_make_unique(join="-")
adata_WP.obs_names_make_unique(join="-")


#%% FROM SEURAT FILES WC

# Read the .mtx file along with the barcodes and gene names
adata_WC = sc.read_mtx("/home/liza/Documents/PhD/scRNAseq_Inbal/WC.mtx")
genes = pd.read_csv("/home/liza/Documents/PhD/scRNAseq_Inbal/WCgenes.tsv", header=None, sep='\t')[0]
barcodes = pd.read_csv("/home/liza/Documents/PhD/scRNAseq_Inbal/WCbarcodes.tsv", header=None, sep='\t')[0]

adata_WC = adata_WC.transpose()

# Check dimensions and assign names
if adata_WC.shape[1] == len(genes):
    # If the number of columns in adata.X matches the number of genes, assign var_names
    adata_WC.var_names = genes
    adata_WC.obs_names = barcodes
else:
    # If there's a mismatch, raise an error or alert
    raise ValueError("Mismatch between the number of genes and the number of columns in the AnnData object.")


# Assign the gene names and cell barcodes
adata_WC.var_names = genes
adata_WC.obs_names = barcodes

# Ensure that var_names and obs_names do not have duplicates
adata_WC.var_names_make_unique(join="-")
adata_WC.obs_names_make_unique(join="-")

#%% If you neas to read bc_matrix files

# adata_GC1 = sc.read_10x_mtx(
#     '/home/liza/Documents/PhD/scRNAseq_Inbal/filtered_feature_bc_matrix - dblGATA_Control', var_names='gene_symbols', cache=True)

# #%%
# # Get gene names from both objects
# genes_adata_GC = adata_GC.var_names
# genes_adata_GC1 = adata_GC1.var_names

# # Find the unique genes in each object and the intersection
# unique_to_adata_GC = set(genes_adata_GC) - set(genes_adata_GC1)
# unique_to_adata_GC1 = set(genes_adata_GC1) - set(genes_adata_GC)
# common_genes = set(genes_adata_GC).intersection(set(genes_adata_GC1))

# # Print the number of unique and common genes
# print(f"Unique to adata_GC: {len(unique_to_adata_GC)}")
# print(f"Unique to adata_GC1: {len(unique_to_adata_GC1)}")
# print(f"Common genes: {len(common_genes)}")

# #%%
# # Get obs names from both objects
# obs_adata_GC = adata_GC.obs_names
# obs_adata_GC1 = adata_GC1.obs_names

# # Find the unique genes in each object and the intersection
# unique_to_adata_GC = set(obs_adata_GC) - set(obs_adata_GC1)
# unique_to_adata_GC1 = set(obs_adata_GC1) - set(obs_adata_GC)
# common_obs = set(obs_adata_GC).intersection(set(obs_adata_GC1))

# # Print the number of unique and common genes
# print(f"Unique to adata_GC: {len(unique_to_adata_GC)}")
# print(f"Unique to adata_GC1: {len(unique_to_adata_GC1)}")
# print(f"Common obs: {len(common_obs)}")


#%%

# adata_GP1 = sc.read_10x_mtx(
#     '/home/liza/Documents/PhD/scRNAseq_Inbal/filtered_feature_bc_matrix - dblGATA_PyMT', var_names='gene_symbols', cache=True)

# adata_WC1 = sc.read_10x_mtx(
#     '/home/liza/Documents/PhD/scRNAseq_Inbal/filtered_feature_bc_matrix - WT_Control', var_names='gene_symbols', cache=True)

# adata_WP1 = sc.read_10x_mtx(
#     '/home/liza/Documents/PhD/scRNAseq_Inbal/filtered_feature_bc_matrix - WT_PyMT', var_names='gene_symbols', cache=True)

# sc.pp.filter_cells(adata_GC1, min_genes=200)
# sc.pp.filter_genes(adata_GC1, min_cells=10)

# sc.pp.filter_cells(adata_GP1, min_genes=200)
# sc.pp.filter_genes(adata_GP1, min_cells=10)

# sc.pp.filter_cells(adata_WC1, min_genes=200)
# sc.pp.filter_genes(adata_WC1, min_cells=10)

# sc.pp.filter_cells(adata_WP1, min_genes=200)
# sc.pp.filter_genes(adata_WP1, min_cells=10)

#%% Put all the data together

adata_GC.obs['sample'] = 'dblGATA_Control'
adata_GP.obs['sample'] = 'dblGATA_PyMT'
adata_WC.obs['sample'] = 'WT_Control'
adata_WP.obs['sample'] = 'WT_PyMT'

adata = sc.concat([adata_GC, adata_GP, adata_WC, adata_WP])

# Subset the adata object for dblGATA samples
adata_dblGATA = adata[adata.obs['sample'].isin(['dblGATA_Control', 'dblGATA_PyMT'])]


# Subset the adata object for WT samples
adata_WT = adata[adata.obs['sample'].isin(['WT_Control', 'WT_PyMT'])]

#%% Check if the values are raw counts

from scipy.sparse import issparse

# Check if the data is stored as a sparse matrix
if issparse(adata.X):
    # Convert to a dense format for viewing
    dense_X = adata.X.toarray()
    print("Adatat is sparse-matrix")
else:
    # If it's already a dense format, just assign it
    dense_X = adata.X

# Now you can view the count values
# For example, view the first 5 genes for the first 5 cells
print(dense_X[:5, :5])

#%% MAKE A COPY of RAW counts

adata_dblGATA.raw = adata_dblGATA
adata_WT.raw = adata_WT
adata.raw = adata

#%% If you need to make sparse matrix from adata

adata.X = csc_matrix(adata.X)  

#%% Filtering the data

adata.var['mt'] = adata.var_names.str.startswith('mt')

adata.var[adata.var.mt == True]
#
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True)

sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')

sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')

upper_lim = np.quantile(adata.obs.n_genes_by_counts.values, .98)
lower_lim = np.quantile(adata.obs.n_genes_by_counts.values, .02)

#

adata = adata[(adata.obs.n_genes_by_counts < upper_lim) & (adata.obs.n_genes_by_counts > lower_lim)]

adata = adata[adata.obs.pct_counts_mt<13, :]

#%% Data normalization

sc.pp.normalize_total(adata, target_sum=1e4)

sc.pp.log1p(adata)

sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

adata = adata[:, adata.var.highly_variable]

sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])

#%% Scaling data

sc.pp.scale(adata, max_value=10)

sc.tl.pca(adata, svd_solver='arpack')

sc.pl.pca(adata, color = 'sample')
#

sc.pl.pca_variance_ratio(adata, log=True)

#%% Plotting UMAP

print(adata.X[:5, :5])

sc.pp.neighbors(adata, n_neighbors=10, n_pcs=30)

sc.tl.umap(adata)

sc.pl.umap(adata)

sc.tl.leiden(adata, resolution = 0.3)
sc.pl.umap(adata, color=['leiden'], frameon = False, legend_loc = 'on data')
sc.pl.umap(adata, color=['sample'], frameon = False)


sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
# sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)


#%% Save TOP ranked genes

top_n_genes = 25

marker_genes = pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(top_n_genes)

marker_genes.to_csv('/home/liza/Documents/PhD/scRNAseq_Inbal/TOPMarkers.csv', sep=',', encoding='utf-8', header='true', index=False)

#%% B-cells - 0, 10

sc.pl.umap(adata, color = ['Cd79a', 'Ms4a1', 'Cd79b', 'Ighd', 'Cd19'], frameon=False, legend_loc = 'on data')

#%% Neutrophils - 1

sc.pl.umap(adata, color = ['G0s2', 'Clec4d', 'S100a9', 'S100a8', 'Cd14'], frameon=False, legend_loc = 'on data')

sc.pl.umap(adata, color=['leiden'], frameon = False, legend_loc = 'on data')

#%% Alveolar macrophage - 6, 8

sc.pl.umap(adata, color = ['Ear2', 'Ear1', 'Cd68', 'Marco', 'Siglecf'], frameon=False, legend_loc = 'on data')

#%% CD8+ T-cell - 3

sc.pl.umap(adata, color = ['Cd8a', 'Cd8b1', 'Il7r', 'Ccr7', 'Klrc1'], frameon=False, legend_loc = 'on data')

#%% Plasmacytoid dendritic cell - 13

sc.pl.umap(adata, color = ['Ms4a6c', 'Plac8', 'Bst2', 'Irf7', 'Irf5'], frameon=False, legend_loc = 'on data')

#%% Dendritic -  5

sc.pl.umap(adata, color = ['Naaa', 'Irf8', 'Cd74', 'Itgax', 'Itgae'], frameon=False, legend_loc = 'on data')

#%% NK Cell - 4

sc.pl.umap(adata, color = ['Nkg7', 'Klra8', 'Klra4', 'Klrb1c'], frameon=False, legend_loc = 'on data')

#%% Interstitial macrophage - 

sc.pl.umap(adata, color = ['C1qc', 'C1qa', 'Pf4', 'Cd74', 'Adgre1'], frameon=False, legend_loc = 'on data')

#%% CD4+ T-cells - 2

sc.pl.umap(adata, color = ['Cd4', 'Tnfrsf4', 'Il7r', 'Ccr7', 'Icos'], frameon=False, legend_loc = 'on data')

#%%  T-reg - 

sc.pl.umap(adata, color = ['Foxp3', 'Ctla4', 'Il2ra', 'Tnfrsf18'], frameon=False, legend_loc = 'on data')

#%%

sc.pl.umap(adata, color = ['Siglecf'], frameon=False, legend_loc = 'on data')


#%% ANNOTATION CELL CLUSTERS

cell_type = {
    '0': "B",
    '1': "Neut",
    '2': "CD4+",
    '3': "CD8+",
    '4': "NK",
    '5': "Dend",
    '6': "Macro",
    '7': "Th1",
    '8': "Macro",
    '9': "Th2",
    '10': "B",
    '11': "Endo",
    '12': "PyMT",
    '13': "pDCs",
    '14': "Pc"
    }

adata.obs['cell_type'] = adata.obs['leiden'].map(cell_type)

sc.pl.umap(adata, color='cell_type', legend_loc = 'on data', title='Cell types', frameon=False)
#%%

import numpy as np
import pandas as pd

# Define a function to process results and return a DataFrame
def process_results(adata_subset):
    results = adata_subset.uns['rank_genes_groups']
    
    out = []
    for group in results['names'].dtype.names:
        for idx in range(len(results['names'][group])):
            gene_name = results['names'][group][idx]
            cell_idx = np.where(adata_subset.var_names == gene_name)[0]  # Get the cell index for the gene_name
            if len(cell_idx) > 0:
                cell_type = adata_subset.obs['cell_type'][cell_idx[0]]  # Extract cell type info using the cell index
                out.append([
                    gene_name,
                    results['scores'][group][idx],
                    results['pvals_adj'][group][idx],
                    results['logfoldchanges'][group][idx],
                    cell_type  # Append the cell_type information here
                ])
            
    markers = pd.DataFrame(out, columns = ['names', 'scores', 'pvals_adj', 'logfoldchanges', 'cell_type'])
    
    markers = markers[(markers.pvals_adj.astype(float) < 0.05) & (abs(markers.logfoldchanges.astype(float)) > 1)]
    
    return markers

# Subset the adata object for dblGATA samples
adata_dblGATA = adata[adata.obs['sample'].isin(['dblGATA_Control', 'dblGATA_PyMT'])]


# Subset the adata object for WT samples
adata_WT = adata[adata.obs['sample'].isin(['WT_Control', 'WT_PyMT'])]

# Process the subsets
markers_dblGATA = process_results(adata_dblGATA)
markers_WT = process_results(adata_WT)

# Write the results to separate CSV files
markers_dblGATA.to_csv('/home/liza/Documents/PhD/scRNAseq_Inbal/Markers_dblGATA.csv', sep=',', encoding='utf-8', header=True, index=False)
markers_WT.to_csv('/home/liza/Documents/PhD/scRNAseq_Inbal/Markers_WT.csv', sep=',', encoding='utf-8', header=True, index=False)

#%% Make count matrix

# Check if it has the todense attribute (i.e., if it's a sparse matrix)
if hasattr(adata_dblGATA.X, 'todense'):
    df1 = pd.DataFrame(adata_dblGATA.X.todense(),
                       index=adata_dblGATA.obs_names,  # Swapped this
                       columns=adata_dblGATA.var_names)  # with this
else:
    df1 = pd.DataFrame(adata_dblGATA.X,
                       index=adata_dblGATA.obs_names,  # Swapped this
                       columns=adata_dblGATA.var_names)  # with this

# Save the DataFrame to a CSV file
df1.T.to_csv('/home/liza/Documents/PhD/scRNAseq_Inbal/count_matrix_dblGATA.tsv', sep='\t', encoding='utf-8', header='true')

if hasattr(adata_WT.X, 'todense'):
    df2 = pd.DataFrame(adata_WT.X.todense(),
                       index=adata_WT.obs_names,  # Swapped this
                       columns=adata_WT.var_names)  # with this
else:
    df2 = pd.DataFrame(adata_WT.X,
                       index=adata_WT.obs_names,  # Swapped this
                       columns=adata_WT.var_names)  # with this


df2.T.to_csv('/home/liza/Documents/PhD/scRNAseq_Inbal/count_matrix_WT.tsv', sep='\t', encoding='utf-8', header='true')

#%% add cell_ID info

adata_dblGATA.obs_names = adata_dblGATA.obs_names.str.replace('-', '.')

# Extract the cell_ID and cell_type into a DataFrame
df_cell_info = pd.DataFrame({
    'cell_ID': adata_dblGATA.obs_names,
    'cell_type': adata_dblGATA.obs['cell_type']
})

# Save to CSV
df_cell_info.to_csv('/home/liza/Documents/PhD/scRNAseq_Inbal/cell_dblGATA_info.tsv', index=False, sep='\t', encoding='utf-8', header='true')

adata_WT.obs_names = adata_WT.obs_names.str.replace('-', '.')

# Extract the cell_ID and cell_type into a DataFrame
df_cell_info = pd.DataFrame({
    'cell_ID': adata_WT.obs_names,
    'cell_type': adata_WT.obs['cell_type']
})

df_cell_info.to_csv('/home/liza/Documents/PhD/scRNAseq_Inbal/cell_WT_info.tsv', index=False, sep='\t', encoding='utf-8', header='true')

#%% Visualize UMAP

sc.pl.umap(adata[adata.obs['sample'].isin(['WT_PyMT', 'dblGATA_PyMT'])], color=['leiden', 'sample'], frameon = False, legend_loc = 'on data')

sc.pl.umap(adata[adata.obs['sample'].isin(['WT_Control', 'dblGATA_Control'])], color=['leiden', 'sample'], frameon = False, legend_loc = 'on data')

#%% CIBERSORTx count_matrix

#Group top genes

sc.tl.rank_genes_groups(adata_dblGATA, groupby='cell_type', method='t-test', n_genes=500)
sc.tl.rank_genes_groups(adata_WT, groupby='cell_type', method='t-test', n_genes=500)

top_genes_dblGATA = {}
for cluster in adata_dblGATA.uns['rank_genes_groups']['names'].dtype.names:
    top_genes_dblGATA[cluster] = adata_dblGATA.uns['rank_genes_groups']['names'][cluster][:500]
    
top_genes_WT = {}
for cluster in adata_WT.uns['rank_genes_groups']['names'].dtype.names:
    top_genes_WT[cluster] = adata_WT.uns['rank_genes_groups']['names'][cluster][:500]
    
#%% Convert the sparse matrix to a dense one if it's in sparse format

import pandas as pd
from scipy.sparse import issparse

if issparse(adata_dblGATA.raw.X):
    raw_counts_dblGATA = adata_dblGATA.raw.X.toarray()
else:
    raw_counts_dblGATA = adata_dblGATA.raw.X

# Create a DataFrame for raw counts
counts_dblGATA = pd.DataFrame(raw_counts_dblGATA, index=adata_dblGATA.raw.obs_names, columns=adata_dblGATA.raw.var_names)

if issparse(adata_WT.raw.X):
    raw_counts_WT = adata_WT.raw.X.toarray()
else:
    raw_counts_WT = adata_WT.raw.X

# Create a DataFrame for raw counts
counts_WT = pd.DataFrame(raw_counts_WT, index=adata_WT.raw.obs_names, columns=adata_WT.raw.var_names)

#%% Add the cell type annotations to the raw counts DataFrame

counts_dblGATA['cell_type'] = adata_dblGATA.obs['cell_type']
counts_WT['cell_type'] = adata_WT.obs['cell_type']

#%% Slice for unique genes

unique_genes_dblGATA = sorted(set([gene for genes in top_genes_dblGATA.values() for gene in genes]))
unique_genes_WT = sorted(set([gene for genes in top_genes_WT.values() for gene in genes]))

# Initialize a DataFrame to hold the averaged raw counts
avg_raw_counts_dblGATA = pd.DataFrame(index=unique_genes_dblGATA)
avg_raw_counts_WT = pd.DataFrame(index=unique_genes_WT)

#%% Loop through each cell type and calculate the mean raw counts for the top genes

for cell_type, genes in top_genes_dblGATA.items():
    cell_type_mask = adata_dblGATA.obs['cell_type'] == cell_type
    
    mean_expression = counts_dblGATA.loc[cell_type_mask, list(genes)].mean(axis=0)
    # Add the mean expression values to the avg_raw_counts_df DataFrame
    avg_raw_counts_dblGATA[cell_type] = mean_expression.reindex(avg_raw_counts_dblGATA.index)

    
for cell_type, genes in top_genes_WT.items():
    cell_type_mask = adata_WT.obs['cell_type'] == cell_type
    
    mean_expression = counts_WT.loc[cell_type_mask, list(genes)].mean(axis=0)
    
    # Add the mean expression values to the avg_raw_counts_df DataFrame
    avg_raw_counts_WT[cell_type] = mean_expression.reindex(avg_raw_counts_WT.index)
    
#%% Get rid of NA and genes with small number of counts

signature_matrix_dblGATA = avg_raw_counts_dblGATA.fillna(0)
signature_matrix_WT = avg_raw_counts_WT.fillna(0)

# Sum the counts across all samples for each gene
gene_counts_sum_dblGATA = signature_matrix_dblGATA.sum(axis=1)
gene_counts_sum_WT = signature_matrix_WT.sum(axis=1)


# Filter genes where the sum of counts is greater than 10
signature_matrix_dblGATA = signature_matrix_dblGATA[gene_counts_sum_dblGATA > 10]
signature_matrix_WT = signature_matrix_WT[gene_counts_sum_WT > 10]


#%% Save the signature matrix to a TSV file

# Get rid of mitochondrial genes:

# Boolean mask to identify mitochondrial genes in dablGATA
mt_genes_mask = signature_matrix_dblGATA.index.str.startswith('mt-')
# Filter out mitochondrial genes
signature_matrix_dblGATA = signature_matrix_dblGATA[~mt_genes_mask]

# Boolean mask to identity genes in WT
mt_genes_mask = signature_matrix_WT.index.str.startswith('mt-')
# Filter out mitochondrial genes
signature_matrix_WT = signature_matrix_WT[~mt_genes_mask]


signature_matrix_dblGATA.to_csv('/home/liza/Documents/PhD/scRNAseq_Inbal/signature_dblGATA.tsv', sep='\t', header=True, index=True)
signature_matrix_WT.to_csv('/home/liza/Documents/PhD/scRNAseq_Inbal/signature_WT.tsv', sep='\t', header=True, index=True)
