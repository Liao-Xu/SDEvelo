# Standard Library Imports
import time

# Third-Party Numerical and Scientific Libraries
import numpy as np
from scipy.stats import kde
import anndata
# from concurrent.futures import ProcessPoolExecutor


# Data Visualization Libraries
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

    
class SimData:
    def __init__(self, total_time=1, dt=0.005, K=10, n_vars=200, seed=0, 
                 a=None, beta=None, c=None, gamma=None, u0=None, s0=None, sig_1=2.0, sig_2=2.0):
        # Simulation Parameters
        self.total_time = total_time
        self.dt = dt
        self.K = K
        self.n_vars = a.shape[0] if a is not None else n_vars
        self.seed = seed

        # Model Parameters
        self.u0 = u0 if u0 is not None else 0.1*np.random.rand(self.n_vars)
        self.s0 = s0 if s0 is not None else 0.1*np.random.rand(self.n_vars)
        # self.u0 = 0.0
        # self.s0 = 0.0
        self.a = a if a is not None else 0.2 + 0.6 *np.random.rand(self.n_vars)
        # self.b = b if b is not None else 100.0
        self.c = c if c is not None else 2 + 10 *np.random.rand(self.n_vars)*np.random.rand(self.n_vars)
        self.beta = beta if beta is not None else 2 + 10 *np.random.rand(self.n_vars)
        self.gamma = gamma if gamma is not None else 2 + 10 *np.random.rand(self.n_vars)

        # Noise Parameters
        self.sigma_1 = sig_1 * np.random.rand(self.n_vars)
        self.sigma_2 = sig_2 * np.random.rand(self.n_vars)
        
        self.adata = None


    def generate(self):
        np.random.seed(self.seed if self.seed is not None else np.random.randint(10000))
        n_obs_per_trajectory = int(self.total_time / self.dt)
        n_obs = n_obs_per_trajectory * self.K

        U_layers = np.zeros((n_obs, self.n_vars))
        S_layers = np.zeros((n_obs, self.n_vars))

        for k in range(self.K):
            U, S = np.zeros((n_obs_per_trajectory, self.n_vars)), np.zeros((n_obs_per_trajectory, self.n_vars))
            Z1, Z2 = np.random.normal(0, 1, self.n_vars), np.random.normal(0, 1, self.n_vars)
            
            U[0] = np.maximum(0 + self.sigma_1 * np.sqrt(self.dt) * Z1, 0)
            S[0] = np.maximum(0 + self.sigma_2 * np.sqrt(self.dt) * Z2, 0)
            for t in range(n_obs_per_trajectory-1):
                time = t * self.dt
                alpha_t = self.c / (1 + np.exp(100.0 * (time - self.a)))
                # alpha_t = self.c / (1 + np.exp(100.0 * (time - 0.5)))
                Z1, Z2 = np.random.normal(0, 1, self.n_vars), np.random.normal(0, 1, self.n_vars)
                
                dU = (alpha_t - self.beta * U[t]) * self.dt + self.sigma_1 * np.sqrt(self.dt) * Z1
                dS = (self.beta * U[t] - self.gamma * S[t]) * self.dt + self.sigma_2 * np.sqrt(self.dt) * Z2
                
                U[t + 1] = np.maximum(U[t] + dU, 0)
                S[t + 1] = np.maximum(S[t] + dS, 0)
            U = self.u0 + U
            S = self.s0 + S

            start, end = k * n_obs_per_trajectory, (k + 1) * n_obs_per_trajectory
            U_layers[start: end], S_layers[start: end] = U, S

        obs = {"true_t": np.tile(np.arange(n_obs_per_trajectory) * self.dt, self.K)}
        var_params = {
            "true_t_": self.a, "true_alpha": self.c, "true_beta": self.beta, 
            "true_gamma": self.gamma, #"true_scaling": self.b, 
            "true_sigma_1": self.sigma_1, "true_sigma_2": self.sigma_2,
            "true_u0": self.u0, "true_s0": self.s0
        }
        layers = {"unspliced": U_layers, "spliced": S_layers}
        
        self.adata = anndata.AnnData(S_layers, obs=obs, var=var_params, layers=layers)
        return self.adata

    def plot_scatter(self, n1, n2, dpi=200):
        spliced = self.adata.layers['spliced']
        unspliced = self.adata.layers['unspliced']
        true_t = self.adata.obs['true_t']

        fig, axes = plt.subplots(n1, n2, figsize=(n2*4, n1*3), dpi=dpi)
        for i in range(n1 * n2):
            ax = axes[i//n2, i%n2] if (n1 > 1 or n2 > 1) else axes
            sc = ax.scatter(spliced[:, i], unspliced[:, i], c=true_t, s=30, alpha=0.5, cmap='viridis_r')
            plt.colorbar(sc, ax=ax)
            # ax.set_title(f'Scatter Plot for Variable {i}')
            ax.set_xlabel('Spliced')
            ax.set_ylabel('Unspliced')

        plt.tight_layout()
        plt.show()