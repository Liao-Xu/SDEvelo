import numpy as np
import scvelo as scv

# Visualization Libraries
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import seaborn as sns


def plot_streamline(adata, args):
    scv.tl.velocity_graph(adata, vkey='sde_velocity', n_jobs=args.scv_n_jobs)
    scv.tl.velocity_confidence(adata, vkey='sde_velocity')
    plt.rcParams.update({'figure.max_open_warning': 0})
    scv.pl.velocity_embedding_stream(adata, 
                                          legend_loc='right_margin', 
                                          vkey="sde_velocity", 
                                          basis=args.vis_key, 
                                          color=[args.vis_type_col], 
                                          title=f"SDEvelo_seed={args.seed}_Epoch={args.nEpochs}", 
                                          cmap="Blues",
                                          figsize=(4, 3), 
                                          dpi=300, 
                                          show=True) 
    
def plot_latent_time(adata, args, figsize=(4.5, 3), dpi=300):
    plt.rcParams.update({'figure.max_open_warning': 0})
    plt.figure(figsize=figsize, dpi=dpi)
    ax = sns.scatterplot(x=adata.obsm[args.vis_key][:, 0], y=adata.obsm[args.vis_key][:, 1],
                         c=adata.obs['latent_time'], cmap="viridis_r")
    norm = plt.Normalize(0, 1)
    sm = plt.cm.ScalarMappable(cmap="viridis_r", norm=norm)
    sm.set_array([])

    # Remove the legend and add a colorbar
    ax.collections[0].set_sizes([150])
    ax.figure.colorbar(sm)
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