import numpy as np
import scvelo as scv
import pandas as pd
import warnings
import scipy.sparse
from matplotlib import MatplotlibDeprecationWarning


# Suppress specific MatplotlibDeprecationWarning from the mentioned module
warnings.filterwarnings("ignore", category=MatplotlibDeprecationWarning, module="scvelo.plotting.utils")

# Visualization Libraries
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import seaborn as sns


def plot_gene_scatter(adata, args, top_genes, nrows=1, ncols=None, figsize=None, dpi=300):
    """
    Plots scatter plots for the specified top genes, showing spliced vs unspliced expression,
    overlaid with predicted values colored by latent time.

    Parameters:
    - adata: AnnData object containing the data.
    - args: Configuration arguments (if needed).
    - top_genes: List of gene names to plot.
    - nrows: Number of rows in the subplot grid.
    - ncols: Number of columns in the subplot grid. If None, it will be calculated based on the number of genes.
    - figsize: Tuple specifying the figure size. If None, it will be set automatically.
    - dpi: Dots per inch (resolution) for the figure. Default is 300.

    Returns:
    - None. Displays the generated plots.
    """
    # Ensure top_genes is a list
    if isinstance(top_genes, str):
        top_genes = [top_genes]
    
    # Determine grid size based on the number of genes
    num_genes = len(top_genes)
    if ncols is None:
        ncols = num_genes if nrows == 1 else (num_genes + nrows - 1) // nrows

    # Set default figure size if not provided
    if figsize is None:
        figsize = (5 * ncols, 5 * nrows)

    # Create subplots with specified dpi
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, dpi=dpi)
    fig.subplots_adjust(hspace=0.4, wspace=0.4)

    # Flatten axes for easy iteration
    axes = np.array(axes).flatten()

    # Fetch the cell type colors from your AnnData object
    celltype_colors = adata.uns.get('celltype_colors', None)
    if celltype_colors is None:
        # Generate a color palette if celltype_colors is not defined
        cell_types = adata.obs['celltype'].unique()
        palette = sns.color_palette('hls', len(cell_types))
        celltype_colors = dict(zip(cell_types, palette))

    for ax, gene in zip(axes, top_genes):
        if gene not in adata.var_names:
            print(f"Warning: Gene '{gene}' not found in adata.var_names. Skipping.")
            continue

        # Extracting the gene's expression data
        gene_idx = adata.var_names.get_loc(gene)
        spliced_layer = adata.layers['Ms']
        unspliced_layer = adata.layers['Mu']
        pred_spliced_layer = adata.layers['pred_s']
        pred_unspliced_layer = adata.layers['pred_u']

        # Handle sparse and dense matrices
        spliced_data = spliced_layer[:, gene_idx].A1 if scipy.sparse.issparse(spliced_layer) else spliced_layer[:, gene_idx]
        unspliced_data = unspliced_layer[:, gene_idx].A1 if scipy.sparse.issparse(unspliced_layer) else unspliced_layer[:, gene_idx]
        pred_spliced_data = pred_spliced_layer[:, gene_idx].A1 if scipy.sparse.issparse(pred_spliced_layer) else pred_spliced_layer[:, gene_idx]
        pred_unspliced_data = pred_unspliced_layer[:, gene_idx].A1 if scipy.sparse.issparse(pred_unspliced_layer) else pred_unspliced_layer[:, gene_idx]
        
        # Extract celltype and latent_time information
        celltype = adata.obs['celltype']
        latent_time = adata.obs['latent_time']

        # Create a DataFrame for easier plotting
        data = pd.DataFrame({
            'Spliced': spliced_data,
            'Unspliced': unspliced_data,
            'Pred_Spliced': pred_spliced_data,
            'Pred_Unspliced': pred_unspliced_data,
            'Celltype': celltype,
            'LatentTime': latent_time
        })

        # Plotting on the current axis with celltype colors
        sns.scatterplot(
            x='Spliced', y='Unspliced', hue='Celltype', data=data,
            palette=celltype_colors, legend=False, edgecolor='none', alpha=0.05, ax=ax
        )

        # Overlay predicted spliced vs. unspliced colored by latent time
        scatter = ax.scatter(
            x=data['Pred_Spliced'], y=data['Pred_Unspliced'],
            c=data['LatentTime'], cmap="viridis_r", edgecolor='none', s=100, alpha=0.7
        )

        # Style adjustments
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.set_title(f'{gene}', fontstyle='italic')
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set(xticklabels=[], yticklabels=[], xlabel=None, ylabel=None)

    # Remove unused subplots if any
    for i in range(num_genes, nrows * ncols):
        fig.delaxes(axes[i])

    # Add colorbar for latent time to the rightmost subplot
    cbar_ax = fig.add_axes([0.92, 0.15, 0.015, 0.7])  # [left, bottom, width, height]
    norm = plt.Normalize(data['LatentTime'].min(), data['LatentTime'].max())
    cbar = plt.colorbar(
        plt.cm.ScalarMappable(norm=norm, cmap='viridis_r'),
        cax=cbar_ax, orientation='vertical', label='Latent Time'
    )

    plt.tight_layout(rect=[0, 0, 0.9, 1])  # Adjust layout to make space for colorbar
    plt.show()

def plot_streamline(adata, args, figsize=(10, 7), dpi=300, **kwargs):
    scv.tl.velocity_graph(adata, vkey='sde_velocity', n_jobs=args.scv_n_jobs)
    scv.tl.velocity_confidence(adata, vkey='sde_velocity')
    plt.rcParams.update({'figure.max_open_warning': 0})
    scv.pl.velocity_embedding_stream(adata, 
                                          legend_loc='right_margin', 
                                          vkey="sde_velocity", 
                                          basis=args.vis_key, 
                                          color=[args.vis_type_col], 
                                          title=f"SDEvelo_seed={args.seed}_Epoch={args.nEpochs}", 
                                          # cmap="Blues",
                                          figsize=figsize, 
                                          dpi=dpi, 
                                          show=True, **kwargs) 
    
def plot_latent_time(adata, args, figsize=(15, 8), dpi=300):
    plt.rcParams.update({'figure.max_open_warning': 0})
    plt.figure(figsize=figsize, dpi=dpi)
    ax = sns.scatterplot(x=adata.obsm[args.vis_key][:, 0], y=adata.obsm[args.vis_key][:, 1],
                         c=adata.obs['latent_time'], cmap="viridis_r")
    norm = plt.Normalize(0, 1)
    sm = plt.cm.ScalarMappable(cmap="viridis_r", norm=norm)
    sm.set_array([])

    # Remove the legend and add a colorbar
    ax.collections[0].set_sizes([150])
    plt.colorbar(sm, ax=ax)
    plt.axis('off')
    plt.show()


def plot_noise_histogram(adata, bins=50, figsize=(6, 3), dpi=300):
    """
    Plot a stacked histogram of sigma1 and sigma2 from the SDEvelo model.

    - param model: The model object containing sigma1 and sigma2 as attributes.
    - param bins: int, the number of bins for the histogram.
    - param figsize: tuple of int, the size of the figure (width, height).
    - param dpi: int, the resolution in dots per inch.
    """
    # Extract sigma1 and sigma2 from the model
    sigma1 = adata.var['fit_sigma_1']
    sigma2 = adata.var['fit_sigma_2']

    # Define the common range for both datasets
    min_bin = min(min(sigma1), min(sigma2))
    max_bin = max(max(sigma1), max(sigma2))
    common_bins = np.linspace(min_bin, max_bin, bins)

    # Create the histogram
    plt.rcParams.update({'figure.max_open_warning': 0})
    plt.figure(figsize=figsize, dpi=dpi)
    plt.hist([sigma1, sigma2], bins=common_bins, stacked=True, color=['blue', 'red'], alpha=0.7)

    plt.xlabel('Value')
    plt.ylabel('Frequency')
    plt.legend(['Sigma1', 'Sigma2'])

    # Save the figure
    # plt.savefig(filename, dpi=dpi, bbox_inches='tight')
    plt.show()
    # plt.close()  # Close the figure to avoid displaying it in notebooks or scripts

def plot_gene_expression(adata, top_genes):
    # Define the grid size
    n1, n2 = 1, len(top_genes)  # Adjust based on how many genes you want to plot

    # Create subplots
    fig, axes = plt.subplots(n1, n2, figsize=(3 * n2, 3))  # Adjust the figure size as needed
    fig.subplots_adjust(hspace=0.4, wspace=0.4)  # Adjust spacing between plots

    # Flatten axes for easy iteration, if there's only one row or column ensure it is iterable
    if n1 * n2 == 1:
        axes = [axes]
    else:
        axes = axes.flatten()

    # Fetch the cell type colors from your AnnData object
    celltype_colors = adata.uns['celltype_colors']

    for ax, gene in zip(axes, top_genes):
        # Extracting the gene's expression data
        spliced_data = adata.layers['Ms'][:, adata.var.index.get_loc(gene)].toarray()
        unspliced_data = adata.layers['Mu'][:, adata.var.index.get_loc(gene)].toarray()
        pred_spliced_data = adata.layers['pred_s'][:, adata.var.index.get_loc(gene)].toarray()
        pred_unspliced_data = adata.layers['pred_u'][:, adata.var.index.get_loc(gene)].toarray()
        
        # Extract celltype and latent_time information
        celltype = adata.obs['celltype']
        latent_time = adata.obs['latent_time']

        # Create a DataFrame for easier plotting
        data = pd.DataFrame({
            'Spliced': spliced_data.ravel(),
            'Unspliced': unspliced_data.ravel(),
            'Pred_Spliced': pred_spliced_data.ravel(),
            'Pred_Unspliced': pred_unspliced_data.ravel(),
            'Celltype': celltype,
            'LatentTime': latent_time
        })

        # Plotting on the current axis with celltype colors
        sns.scatterplot(x='Spliced', y='Unspliced', hue='Celltype', data=data, palette=celltype_colors, legend=False, edgecolor='none', alpha=0.05, ax=ax)

        # Overlay predicted spliced vs. unspliced colored by latent time
        scatter = ax.scatter(x=data['Pred_Spliced'], y=data['Pred_Unspliced'], c=data['LatentTime'], cmap="viridis_r", edgecolor='none', s=100, alpha=0.7)

        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.set_title(f'{gene}', fontstyle='italic')
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set(xticklabels=[], yticklabels=[], xlabel=None, ylabel=None)

    # If you end up with an extra subplot, you can hide it
    if len(top_genes) < n1 * n2:
        for i in range(len(top_genes), n1 * n2):
            fig.delaxes(axes[i])
    plt.tight_layout()
    plt.show()