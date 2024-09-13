import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from statsmodels.stats.multitest import multipletests

def infer_gene_correlations(adata, significance_threshold=0.05):
    """
    Analyze correlations between gene expressions and latent time in an AnnData object.

    Parameters:
    adata (AnnData): An AnnData object containing latent time and gene expression data.
    significance_threshold (float): The threshold for statistical significance.

    Returns:
    AnnData: The input AnnData object with additional annotations for gene correlations.
    """
    # Ensure that the latent time is a numpy array and handle potential missing values
    latent_time = adata.obs['latent_time'].dropna().values

    # Handle potential missing values in gene expression data
    expression = adata.X.A
    valid_indices = ~np.isnan(expression).any(axis=1) & ~np.isnan(latent_time)
    expression = expression[valid_indices]
    latent_time = latent_time[valid_indices]

    # Vectorized calculation of correlations and p-values
    correlations, p_values = zip(*[pearsonr(expression[:, j], latent_time) for j in range(expression.shape[1])])

    # Adding the results to the AnnData object
    adata.var['correlation'] = pd.Series(correlations, index=adata.var_names)
    adata.var['correlation_abs'] = pd.Series(np.abs(correlations), index=adata.var_names)
    adata.var['p_value'] = pd.Series(p_values, index=adata.var_names)

    # Applying Benjamini-Hochberg procedure
    p_adjusted = multipletests(p_values, method='fdr_bh')[1]
    adata.var['p_value_adj'] = pd.Series(p_adjusted, index=adata.var_names)

    # Sorting genes and dropping unnecessary columns
    genes_sorted = adata.var.sort_values(by='correlation_abs', ascending=False)

    # Filter out significant genes
    significant_genes = genes_sorted[genes_sorted['p_value_adj'] <= significance_threshold]
    
    return adata, significant_genes