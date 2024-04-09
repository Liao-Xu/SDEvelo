import torch
import torch.nn as nn
import numpy as np
import scvelo as scv

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
        self.ukey = args.ukey
        self.skey = args.skey
        self.mode = args.mode
        self.device = self.setup_environment(args)
        if self.ukey == 'Mu':
            # Filter and normalize
            scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=args.n_gene)
            # Compute moments
            scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
            
        self.adata = adata
        paras = self.init_para(self.adata)

        if paras.gamma_init is None:
            paras.gamma_init = torch.ones(K)
            
        if paras.a_init is None:
            paras.a_init = torch.ones(K)
            
        if paras.c_init is None:
            paras.c_init = torch.ones(K)
            
        self.K = paras.K
        self.a = torch.tensor(paras.a_init).to(self.device).requires_grad_(True) #original time
        self.b = (100.0*torch.ones(paras.K)).to(self.device).requires_grad_(True)
        self.c = torch.tensor(paras.c_init).to(self.device).requires_grad_(True) #original alpha
        self.h = 0.01
        self.sigma1 = (0.01*torch.ones(paras.K)).to(self.device).requires_grad_(True)
        self.sigma2 = (0.01*torch.ones(paras.K)).to(self.device).requires_grad_(True)
        self.u = Tensor(paras.u_train, paras.K)
        self.s = Tensor(paras.s_train, paras.K)

        train_dat = TensorDataset(self.s, self.u)
        self.trainLoader = DataLoader(train_dat, batch_size=200, shuffle=True)

        
        
        # settings
        if paras.s0_init is None:
            self.s0 = torch.zeros(paras.K).to(self.device)
        if paras.u0_init is None:
            self.u0 = torch.zeros(paras.K).to(self.device)
            
        self.s0 = torch.tensor(paras.s0_init).to(self.device)
        self.u0 = torch.tensor(paras.u0_init).to(self.device)
        self.s_shift = torch.zeros(paras.K).to(self.device).requires_grad_(True)
        self.u_shift = torch.zeros(paras.K).to(self.device).requires_grad_(True)

        self.beta = torch.tensor(paras.c_init).to(self.device).requires_grad_(True)
        self.gamma = torch.tensor(paras.c_init).to(self.device).requires_grad_(True)

        self.optimizer_1 = torch.optim.Adam([self.c, self.beta], 0.5)
        self.optimizer_2 = torch.optim.Adam([self.gamma], 0.1)
        self.optimizer_3 = torch.optim.Adam([self.a, self.sigma1, self.sigma2, self.u_shift, self.s_shift], 1e-3)
        
        self.scheduler_1 = MultiStepLR(self.optimizer_1, milestones=[200], gamma=0.1)
        self.scheduler_2 = MultiStepLR(self.optimizer_2, milestones=[200], gamma=0.1)

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
        Generates SDE paths for u and s variables using Euler-Maruyama method.
        """
        a = self.a
        b = self.b
        c = self.c
        beta = self.beta
        gamma = self.gamma
        h = self.h
        sig1 = self.sigma1*np.sqrt(h)
        sig2 = self.sigma2*np.sqrt(h)
        t_gen = torch.arange(0, 1, h).to(self.device)

        u_sde = torch.zeros(len(t_gen), self.K).to(self.device)
        s_sde = torch.zeros(len(t_gen), self.K).to(self.device)

        u_sde[0] = torch.relu(self.u0 + self.u_shift)
        s_sde[0] = torch.relu(self.s0 + self.s_shift)
        
        for i in range(0, len(t_gen) - 1):
            exp_term = torch.clamp(b * (t_gen[i] - a), min=-50, max=50)
            u_sde_temp = u_sde[i] + h * (c / (1 + torch.exp(exp_term) + 1e-9) - beta * u_sde[i].clone()) + sig1 * torch.randn(self.K).to(self.device)
            s_sde_temp = s_sde[i] + h * (beta * u_sde[i].clone() - gamma * s_sde[i].clone()) + sig2 * torch.randn(self.K).to(self.device)

            u_sde[i + 1] = torch.relu(u_sde_temp)
            s_sde[i + 1] = torch.relu(s_sde_temp)

        return u_sde , s_sde 

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
            if (epoch) % 20 == 0:
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
        Calculate latent time.
        """
        s_raw = self.adata.layers[self.skey]
        u_raw = self.adata.layers[self.ukey]
        u_pred, s_pred = self.generate()
        pred_data = np.concatenate((u_pred, s_pred), axis=1)
        raw_data = np.concatenate((u_raw, s_raw), axis=1)
    
        distances = cdist(raw_data, pred_data, 'sqeuclidean')
        nearest_indices = np.argmin(distances, axis=1)
        t_pred = np.linspace(0, 1, int(1 / self.h))
        latent_time = t_pred[nearest_indices]
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

