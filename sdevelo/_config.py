class Config:
    def __init__(self):
        # Size of mini-batches for training.
        self.batchSz = 200
        
        # Identifier for the CUDA device to use. Default is 0.
        self.cuda_device = 2
        
        # The number of genes to consider in the model. Default is 2000.
        self.n_gene = 2000

        # The number of genes to consider in the model. Default is 2000.
        self.sde_mode = 'original'
        
        # Key for accessing visualization data. Default is 'X_pca'.
        self.vis_key = 'X_pca'
        
        # Determines the coloring scheme for visualization. Default is 'Clusters'.
        self.vis_type_col = 'true_t'

        self.ukey = 'Mu'
        
        self.skey = 'Ms'
        
        # Flag to disable CUDA even if it's available. Default is False (CUDA enabled if available).
        self.no_cuda = False

        self.time_mode = 0 # Default is 0, nearest neigbour for inferring latent time
        
        # Number of parallel jobs for scVelo (if applicable). Default is 10.
        self.scv_n_jobs = 10
        
        # The number of cells to include in the model. Default is 2000.
        self.n_cell = 5000

        self.process = True
        
        # Seed for random number generators to ensure reproducibility. Default is 0.
        self.seed = 0
        
        # The number of training epochs. Default is 50.
        self.nEpochs = 300

        # The mode to calculate velocity.
        self.mode = 2


    