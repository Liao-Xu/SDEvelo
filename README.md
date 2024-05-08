# SDEvelo: a deep generative approach for transcriptional dynamics with cell-specific latent time and multivariate stochastic modeling

## Overview

`sdevelo` leverages advanced stochastic differential equations (SDE) to provide a novel approach to RNA velocity analysis in single-cell RNA sequencing (scRNA-seq). This deep generative model accurately captures the complex, stochastic nature of transcriptional dynamics, offering new insights into cell differentiation and state transitions.


## System Requirements

- **Operating Systems**: Linux (Ubuntu, CentOS), macOS, Windows 10.
- **Python Version**: Python 3.6 and above.
- **Dependencies**:
  anndata==0.10.7
  matplotlib==3.7.1
  numpy==1.23.5
  scipy==1.8.1
  scvelo==0.2.5
  seaborn==0.11.2
  torch==1.13.1+cu117
- **Hardware Requirements**: No non-standard hardware required.
- **Installation Time**: `sdevelo`'s installation should be completed within approximately 5 minutes.

## Installation Guide

- **Step 1**: Ensure Python 3.6+ is installed on your system.
- **Step 2**: Install `sdevelo`  via pip:

  ```bash
  pip install sdevelo
  ```


## Demo

Experience the power of SDEvelo through our interactive demo provided as a Jupyter Notebook.

### Running the Demo

1. Navigate to the `docs/demo_simulation` directory within this repository.
2. Locate the Jupyter Notebook titled `demo_simulation.ipynb`.
3. Open the notebook in Jupyter Lab or Jupyter Notebook and execute the cells in order.

### Expected Output

By running the demo, you will generate:

- A streamline plot depicting the transcriptional dynamics.
- A latent time heatmap that visualizes the progression of cells over time.

### Expected Run Time

On a typical desktop computer, the demo should complete within approximately 300 seconds.

### Future Demos

We continuously strive to enhance SDEvelo. Stay tuned for additional demos by checking our repository for updates.


## Instruction to use

- **Step 1**: Configure the arguments and parameters for your dataset. Refer to the provided `demo_simulation.ipynb` for examples of data configuration, model execution, and visualization.
- **Step 2**: Run the SDEvelo model.
- **Step 3**: Visualize the results based on the estimated SDEvelo model.
