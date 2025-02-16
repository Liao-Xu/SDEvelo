import torch
import torch.nn as nn
import numpy as np
import scvelo as scv
import scanpy as sc
import importlib
import subprocess
import sys
from packaging import version
from dataclasses import dataclass
from collections import OrderedDict
from torch.utils.data import DataLoader, TensorDataset
from torch.optim.lr_scheduler import MultiStepLR
from scvelo.core import LinearRegression
from scipy.spatial.distance import cdist

@dataclass
class SDEPara:
    K: int
    u_train: np.ndarray
    s_train: np.ndarray
    gamma_init: np.ndarray
    a_init: np.ndarray
    c_init: np.ndarray
    s0_init: np.ndarray
    u0_init: np.ndarray
    
def guassian_kernel(self, source, target, kernel_mul=2.0, kernel_num=5, fix_sigma=None):
    """
    Calculate the Gram kernel matrix.
    - source: A data matrix with dimensions (sample_size_1 x feature_size), representing the first sample set.
    - target: A data matrix with dimensions (sample_size_2 x feature_size), representing the second sample set.
    - kernel_mul: This concept is somewhat unclear but appears to be related to calculating the bandwidth for each kernel.
    - kernel_num: The number of kernels, indicating the use of multiple kernels.
    - fix_sigma: Indicates whether to use a fixed standard deviation.
    
    Returns: A matrix of size ((sample_size_1 + sample_size_2) x (sample_size_1 + sample_size_2)), 
             structured as:
                            [   K_ss K_st
                                K_ts K_tt ]
             where K_ss and K_tt are the kernel matrices within the same sample sets, and K_st and K_ts are the 
             kernel matrices between the two different sample sets.
    """

    n_samples = int(source.size()[0])+int(target.size()[0])
    total = torch.cat([source, target], dim=0) # merge

    total0 = total.unsqueeze(0).expand(int(total.size(0)), \
                                       int(total.size(0)), \
                                       int(total.size(1))).to(self.device)
    total1 = total.unsqueeze(1).expand(int(total.size(0)), \
                                       int(total.size(0)), \
                                       int(total.size(1))).to(self.device)
    L2_distance = ((total0-total1)**2).sum(2) 

    # Calculate the bandwidth for each kernel in a multi-kernel setup.
    if fix_sigma:
        bandwidth = fix_sigma
    else:
        bandwidth = torch.sum(L2_distance.data) / (n_samples**2-n_samples)
        
    bandwidth /= kernel_mul ** (kernel_num // 2)
    bandwidth_list = [bandwidth * (kernel_mul**i) for i in range(kernel_num)]

    # The formula for the Gaussian kernel: exp(-|x-y|/bandwidth)
    kernel_val = [torch.exp(-L2_distance / bandwidth_temp) for \
                  bandwidth_temp in bandwidth_list]

    return sum(kernel_val) # Combine multiple kernels together.


def mmd(self, source, target,  kernel_mul=2.0, kernel_num=5, fix_sigma=None):
    """
    Calculates the Maximum Mean Discrepancy (MMD) between 'source' and 'target' distributions.
    
    Parameters:
    - source: Tensor (sample_size_1 x feature_size) from the source distribution.
    - target: Tensor (sample_size_2 x feature_size) from the target distribution.
    - kernel_mul: Bandwidth multiplier for Gaussian kernels (default=2.0).
    - kernel_num: Number of Gaussian kernels (default=5).
    - fix_sigma: Fixed standard deviation for all kernels; if None, adaptive (default=None).
    
    This function computes the MMD using Gaussian kernels, providing a measure of distance between the distributions of 'source' and 'target'. It calculates kernel matrices for within and between the distributions, normalizes these matrices, and combines them to form the MMD loss.
    
    Returns:
    - The square root of the MMD loss, a scalar indicating the distributional distance.
    """

    n = int(source.size()[0])
    m = int(target.size()[0])

    kernels = guassian_kernel(self, source, target,
                              kernel_mul=kernel_mul, kernel_num=kernel_num, fix_sigma=fix_sigma)
    
    XX = kernels[:n, :n] 
    YY = kernels[n:, n:]
    XY = kernels[:n, n:]
    YX = kernels[n:, :n]
    
    XX = torch.div(XX, n * n).sum(dim=1).view(1,-1)  # K_ss matrix，Source<->Source
    XY = torch.div(XY, -n * m).sum(dim=1).view(1,-1) # K_st matrix，Source<->Target

    YX = torch.div(YX, -m * n).sum(dim=1).view(1,-1) # K_ts matrix,Target<->Source
    YY = torch.div(YY, m * m).sum(dim=1).view(1,-1)  # K_tt matrix,Target<->Target
    
    loss = (XX + XY).sum() + (YX + YY).sum()
    return loss.sqrt()

def Tensor(data, k = 1):
    return torch.tensor(data).view(-1, k).float()

class SDENN:
    """
    Implements the Stochastic Differential Equation Neural Network (SDENN) for modeling and simulation of dynamical systems influenced by stochastic processes. This physics-guided neural network integrates physical laws into the learning process, enhancing prediction accuracy and interpretability.
    
    Parameters:
    - paras: A configuration object containing initialization parameters for the model, including:
        - gamma_init: Initial value(s) for the gamma parameter(s), defaulting to ones if not provided.
        - a_init: Initial value(s) for the 'a' parameter(s), representing original time, defaulting to ones if not provided.
        - c_init: Initial value(s) for the 'c' parameter(s), representing original alpha, defaulting to ones if not provided.
        - K: The number of different processes or variables to model.
        - device: The computational device (CPU, GPU) to use.
        - u_train, s_train: Training data for the 'u' and 's' variables.
        - s0_init, u0_init: Initial conditions for 's' and 'u' variables.
        Additional parameters include learning rates and batch size for the training data loader.
    
    The SDENN class encapsulates the functionality for generating data paths based on Stochastic Differential Equations (SDEs) using the Euler-Maruyama method, training the model to fit observed data, and generating predictions. It employs a combination of adjustable parameters and optimization techniques to closely model the underlying dynamics of the system being studied.
    
    Methods:
    - sde_gen: Generates SDE paths for 'u' and 's' variables.
    - train: Trains the model over a specified number of iterations to minimize the discrepancy between generated and observed data.
    - generate: Produces and returns detached SDE paths for 'u' and 's' variables, based on the trained model parameters.
    
    This model is particularly useful for systems where the underlying dynamics are partially known and can be described by differential equations with stochastic components.
    """

    def __init__(self, args, adata):
        # Initialize model parameters from 'paras' with defaults where applicable.
        self.device = self.setup_environment(args)
        self.ukey = args.ukey
        self.skey = args.skey
        self.mode = args.mode
        self.time_mode = args.time_mode
        self.sde_mode = args.sde_mode
        self.batchSz = args.batchSz
        if args.process == True:
            # Filter and normalize
            scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=args.n_gene)
            scvelo_version = version.parse(scv.__version__)

            # Now we can use the version to determine which code to run
            if scvelo_version < version.parse("0.3.0"):
                print(f"Using scVelo version {scv.__version__}")
                scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
            else:
                print(f"Using scVelo version {scv.__version__}")
                sc.pp.pca(adata)
                sc.pp.neighbors(adata, n_pcs=30, n_neighbors=30)
                scv.pp.moments(adata, n_pcs=None, n_neighbors=None)

        
            
        self.adata = adata
        
        if len(self.adata) > args.n_cell:
            indexes = self.process_data(args)
            paras = self.init_para(self.adata[indexes])
        else:
            paras = self.init_para(self.adata)

        if paras.gamma_init is None:
            paras.gamma_init = torch.ones(K)
            
        if paras.a_init is None:
            paras.a_init = torch.ones(K)
            
        if paras.c_init is None:
            paras.c_init = torch.ones(K)
            
        self.K = paras.K
        self.a = torch.tensor(paras.a_init).to(self.device).requires_grad_(True) #original time
        self.b = (100.0*torch.ones(self.K)).to(self.device).requires_grad_(True)
        self.c = torch.tensor(paras.c_init).to(self.device).requires_grad_(True) #original alpha
        self.h = 0.01
        self.sigma1 = (0.01*torch.ones(self.K)).to(self.device).requires_grad_(True)
        self.sigma2 = (0.01*torch.ones(self.K)).to(self.device).requires_grad_(True)
        self.u = Tensor(paras.u_train, self.K)
        self.s = Tensor(paras.s_train, self.K)

        train_dat = TensorDataset(self.s, self.u)
        self.trainLoader = DataLoader(train_dat, batch_size=self.batchSz, shuffle=True)
        print(len(train_dat))

        
        
        # settings
        if paras.s0_init is None:
            self.s0 = torch.zeros(self.K).to(self.device)
        if paras.u0_init is None:
            self.u0 = torch.zeros(self.K).to(self.device)
            
        self.s0 = torch.tensor(paras.s0_init).to(self.device)
        self.u0 = torch.tensor(paras.u0_init).to(self.device)
        self.s_shift = torch.zeros(self.K).to(self.device).requires_grad_(True)
        self.u_shift = torch.zeros(self.K).to(self.device).requires_grad_(True)

        self.beta = torch.tensor(paras.c_init).to(self.device).requires_grad_(True)
        self.gamma = torch.tensor(paras.c_init).to(self.device).requires_grad_(True)

        self.optimizer_1 = torch.optim.Adam([self.c, self.beta], 0.5)
        self.optimizer_2 = torch.optim.Adam([self.gamma], 0.1)
        self.optimizer_3 = torch.optim.Adam([self.a, self.sigma1, self.sigma2, self.u_shift, self.s_shift], 1e-3)
        
        self.scheduler_1 = MultiStepLR(self.optimizer_1, milestones=[200], gamma=0.1)
        self.scheduler_2 = MultiStepLR(self.optimizer_2, milestones=[200], gamma=0.1)

    def process_data(self, args):
        """
        Process the given adata object using Scanpy and scVelo.
    
        Parameters:
        adata (AnnData): Annotated data matrix.
        args (argparse.Namespace): Arguments, should have 'n_raw_gene' attribute.
        
        Returns:
        indexes (numpy.ndarray): Sorted random indices used for subsetting adata.
        """
    
        # Generate random indices
        num_rows = args.n_cell
        total_rows = len(self.adata)
        random_indices = np.random.choice(total_rows, num_rows, replace=False)
        indexes = np.sort(random_indices)
        return indexes
        
    def setup_environment(self, args):
        """
        Configure the environment for PyTorch and NumPy based on the provided arguments.
        
        Parameters:
        args (argparse.Namespace): The arguments to use for setting up the environment.
                                   It should have the following attributes:
                                   - seed: The random seed.
                                   - no_cuda: Boolean flag to disable CUDA.
                                   - cuda_device: The CUDA device ID to use.
                                   
        Returns:
        torch.device: The device to be used for computation (either CUDA or CPU).
        """
        # Set random seeds for reproducibility
        torch.manual_seed(args.seed)
        np.random.seed(args.seed)
    
        # Determine if CUDA is available and should be used
        args.cuda = not args.no_cuda and torch.cuda.is_available()
        device = torch.device("cuda" if args.cuda else "cpu")
        print(device)
        # Additional configurations if CUDA is being used
        if args.cuda:
            torch.cuda.manual_seed(args.seed)
            torch.cuda.set_device(args.cuda_device)
            torch.backends.cudnn.deterministic = True
    
        return device

    def init_para(self, adata):
        """
        Initialize the SDE parameters based on given s_data and u_data.
    
        Parameters:
        s_data (np.ndarray): s data array.
        u_data (np.ndarray): u data array.
        device (torch.device): Torch device (e.g., 'cpu', 'cuda').
    
        Returns:
        SDEPara: A dataclass object containing initialized SDE parameters.
        """
    
        u_data = adata.layers[self.ukey]
        s_data = adata.layers[self.skey]
        n_cell, n_gene = s_data.shape
        
        # Linear regression
        lr = LinearRegression(fit_intercept=True, percentile=[5, 95])
        lr.fit(s_data, u_data)
    
        # Initial predictions and parameters
        offset_pred = lr.intercept_
        gamma_pred = lr.coef_
        a_init = np.sum(u_data > gamma_pred * s_data, 0) / n_cell
        c_init = np.max(u_data, 0)
    
        # Compute distances and initial conditions
        distances = np.linalg.norm(np.stack((s_data, u_data)), axis=0)
        min_distance_idx = np.argmin(distances, axis=0)
        s0_init = np.array([s_data[min_distance_idx[j], j] for j in range(n_gene)])
        u0_init = np.array([u_data[min_distance_idx[j], j] for j in range(n_gene)])
    
        # Create SDEPara object
        paras = SDEPara(
            K = n_gene,
            s_train = s_data,
            u_train = u_data,
            gamma_init = gamma_pred,
            a_init = a_init,
            c_init = c_init,
            s0_init = s0_init,
            u0_init = u0_init
        )
    
        return paras


    def sde_gen(self):
        """
        Generates SDE paths for u and s variables using either the original Euler-Maruyama method
        or the torchsde method.
    
        Args:
            method (str): The method to use for SDE generation. 
                          Options are 'original' (default) or 'torchsde'.
    
        Returns:
            tuple: (u_sde, s_sde) - The generated SDE paths for u and s variables.
        """
        if self.sde_mode == 'original':
            return self.sde_gen_original()
        elif self.sde_mode == 'torchsde':
            # Check if torchsde is installed
            if importlib.util.find_spec("torchsde") is None:
                print("torchsde is not installed. Attempting to install...")
                try:
                    subprocess.check_call([sys.executable, "-m", "pip", "install", "torchsde"])
                    print("torchsde has been successfully installed.")
                except subprocess.CalledProcessError:
                    print("Failed to install torchsde. Please install it manually.")
                    return self.sde_gen_original()  # Fallback to original method
    
            # If installation was successful or package was already installed, import and use it
            try:
                import torchsde
                return self.sde_gen_torchsde()
            except ImportError:
                print("Failed to import torchsde. Falling back to original method.")
                return self.sde_gen_original()
        else:
            raise ValueError("Invalid method. Choose either 'original' or 'torchsde'.")
    
    def sde_gen_original(self):
        """
        Generates SDE paths for u and s variables using Euler-Maruyama method.
        """
        a = self.a
        b = self.b
        c = self.c
        beta = self.beta
        gamma = self.gamma
        h = self.h
        sig1 = self.sigma1 * np.sqrt(h)
        sig2 = self.sigma2 * np.sqrt(h)
        t_gen = torch.arange(0, 1, h, device=self.device)
        u_sde = torch.zeros(len(t_gen), self.K, device=self.device)
        s_sde = torch.zeros(len(t_gen), self.K, device=self.device)
        u_sde[0] = torch.relu(self.u0 + self.u_shift)
        s_sde[0] = torch.relu(self.s0 + self.s_shift)
        for i in range(0, len(t_gen) - 1):
            exp_term = torch.clamp(b * (t_gen[i] - a), min=-50, max=50)
            denom = 1 + torch.exp(exp_term) + 1e-9
            u_sde_temp = (
                u_sde[i]
                + h * (c / denom - beta * u_sde[i].clone())
                + sig1 * torch.randn(self.K, device=self.device)
            )
            s_sde_temp = (
                s_sde[i]
                + h * (beta * u_sde[i].clone() - gamma * s_sde[i].clone())
                + sig2 * torch.randn(self.K, device=self.device)
            )
            u_sde[i + 1] = torch.relu(u_sde_temp)
            s_sde[i + 1] = torch.relu(s_sde_temp)
        return u_sde, s_sde
    
    def sde_gen_torchsde(self):
        """
        Generates SDE paths for u and s variables using torchsde's sdeint method.
        """
        import torchsde
        class SDEFunc(torchsde.SDEIto):
            def __init__(self, a, b, c, beta, gamma, sigma1, sigma2):
                super().__init__(noise_type='diagonal')
                self.a = a
                self.b = b
                self.c = c
                self.beta = beta
                self.gamma = gamma
                self.sigma1 = sigma1
                self.sigma2 = sigma2
    
            def f(self, t, x):
                u = x[:, 0]
                s = x[:, 1]
                t = t.to(x.device)
                exp_term = torch.clamp(self.b * (t - self.a), min=-50, max=50)
                denom = 1 + torch.exp(exp_term) + 1e-9
                du_dt = self.c / denom - self.beta * u
                ds_dt = self.beta * u - self.gamma * s
                return torch.stack([du_dt, ds_dt], dim=1)
    
            def g(self, t, x):
                batch_size = x.size(0)
                device = x.device
                g_matrix = torch.zeros(batch_size, 2, device=device)
                g_matrix[:, 0] = self.sigma1
                g_matrix[:, 1] = self.sigma2
                return g_matrix
    
        t_span = torch.arange(0, 1, self.h, device=self.device)
        u0 = torch.relu(self.u0 + self.u_shift)
        s0 = torch.relu(self.s0 + self.s_shift)
        x0 = torch.stack([u0, s0], dim=1)
    
        sde_func = SDEFunc(self.a, self.b, self.c, self.beta, self.gamma, self.sigma1, self.sigma2).to(self.device)
    
        sol = torchsde.sdeint(sde_func, x0, t_span, method='euler', dt=self.h)
    
        u_sde = torch.relu(sol[:, :, 0].transpose(0, 1))
        s_sde = torch.relu(sol[:, :, 1].transpose(0, 1))
        return u_sde.T, s_sde.T

    def core(self):
        """
        Trains the model over a specified number of iterations (nIter).
        """
        for epoch in range(self.nIter):
            for _, (dat_s, dat_u) in enumerate(self.trainLoader):
                gen_u, gen_s = self.sde_gen()
                gen_data = torch.cat((gen_u, gen_s), dim=1)
                real_data = torch.cat((dat_u, dat_s), dim=1).to(self.device)
                
                # Compute loss between generated and real data.
                loss = mmd(self, gen_data, real_data)

                # Zero gradients, perform backward pass, and update weights.
                self.optimizer_1.zero_grad()
                self.optimizer_2.zero_grad()
                self.optimizer_3.zero_grad()
                loss.backward()
                self.optimizer_1.step()
                self.optimizer_2.step()
                self.optimizer_3.step()

            # Adjust learning rate.
            self.scheduler_1.step()
            self.scheduler_2.step()

            # Parameter constraints enforcement.
            with torch.no_grad():
                self.a.clamp_(0.1, 1)
                self.u_shift.clamp_(-0.1, 0.1)
                self.s_shift.clamp_(-0.1, 0.1)
                self.c.clamp_(0.01, 1000)
                self.beta.clamp_(0.01, 1000)
                self.gamma.clamp_(0.01, 100)
                self.sigma1.clamp_(0.01, 10.0)
                self.sigma2.clamp_(0.01, 10.0)

            # Optional: Log training progress.
            if (epoch) % 50 == 0:
                # print(epoch+1)
                print(
                'Epoch: %d, Loss: %.3f, alpha: %.2f, beta: %.2f, gamma: %.2f, s1: %.3f, s2: %.3f,  t_m: %.3f,  u_shift: %.3f,  s_shift: %.3f' % 
                (
                    epoch, 
                    loss.item(), 
                    self.c[0].item(), 
                    self.beta[0].item(), 
                    self.gamma[0].item(),
                    self.sigma1[0].item(),
                    self.sigma2[0].item(),
                    self.a[0].item(),
                    self.u_shift[0].item(),
                    self.s_shift[0].item()
                )
            )

    def train(self, nIter):
        """
        Trains the model over a specified number of iterations (nIter).
        """
        self.nIter = nIter
        self.core()
        self.adata = self.velocity_cal()
        self.adata = self.pass_var()
        self.adata = self.latent_time()
        return self.adata
        
    def latent_time(self):
        """
        Calculate latent time for each cell in the dataset.
    
        Uses either Euclidean distance (time_mode=0) or optimal transport (time_mode=1)
        to find the nearest predicted timepoint for each cell.
    
        Returns:
            anndata.AnnData: Updated with 'pred_u', 'pred_s' layers and 'latent_time' observation.
    
        Raises:
            ImportError: If time_mode=1 and 'ot' package is not installed.
            ValueError: If invalid time_mode is provided.
        """
        s_raw = self.adata.layers[self.skey]
        u_raw = self.adata.layers[self.ukey]
        u_pred, s_pred = self.generate()
        pred_data = np.concatenate((u_pred, s_pred), axis=1)
        raw_data = np.concatenate((u_raw, s_raw), axis=1)

        if self.time_mode==0:
            distances = cdist(raw_data, pred_data, 'sqeuclidean')
            nearest_indices = np.argmin(distances, axis=1)
        elif self.time_mode == 1:
            try:
                import ot
            except ImportError:
                raise ImportError("Please install the 'ot' package for optimal transport. You can install it using 'pip install POT'.")

            # Compute the squared Euclidean distance matrix for optimal transport
            cost_matrix = ot.dist(raw_data, pred_data, metric='sqeuclidean')
    
            # Compute the optimal transport plan using the Sinkhorn algorithm
            epsilon = 0.5  # Regularization parameter
            transport_plan = ot.sinkhorn(np.ones(raw_data.shape[0]) / raw_data.shape[0],
                                         np.ones(pred_data.shape[0]) / pred_data.shape[0],
                                         cost_matrix, epsilon)
            
            # Compute the optimal assignment from the transport plan
            nearest_indices = np.argmax(transport_plan, axis=1)
            
        t_pred = np.linspace(0, 1, int(1 / self.h))
        latent_time = t_pred[nearest_indices]
        self.adata.layers['pred_u'] = u_pred[nearest_indices]
        self.adata.layers['pred_s'] = s_pred[nearest_indices]
        self.adata.obs['latent_time'] = latent_time
        return self.adata
        
    
    def pass_var(self):
        self.adata.var['fit_alpha'] = self.c.cpu().detach().numpy()
        self.adata.var['fit_beta'] = self.beta.cpu().detach().numpy()
        self.adata.var['fit_gamma'] = self.gamma.cpu().detach().numpy()
        self.adata.var['fit_t_'] = self.a.cpu().detach().numpy()
        self.adata.var['fit_sigma_1'] = self.sigma1.cpu().detach().numpy()
        self.adata.var['fit_sigma_2'] = self.sigma2.cpu().detach().numpy()
        return self.adata

    def velocity_cal(self):
        beta = self.beta.cpu().detach().numpy()
        gamma = self.gamma.cpu().detach().numpy()
        sigma = self.sigma2.cpu().detach().numpy()
        noise = 0
        h = self.h
        if self.mode==0:
            u = self.adata.layers[self.ukey]
            s = self.adata.layers[self.skey]
        elif self.mode==1:
            u, s = self.generate()
        elif self.mode==2:
            u = self.adata.layers[self.ukey]
            s = self.adata.layers[self.skey]
            noise = sigma*np.sqrt(h)*np.random.randn(beta.shape[0])
        elif self.mode==3:
            u, s = self.generate()
            noise = sigma*np.sqrt(h)*np.random.randn(beta.shape[0])
        self.adata.layers['sde_velocity'] = (beta*u - gamma*s)*h + noise
        return self.adata

    def generate(self):
        """
        Generates and returns detached SDE paths for u and s variables.
        """
        u, s = self.sde_gen()
        return u.detach().cpu().numpy(), s.detach().cpu().numpy()

