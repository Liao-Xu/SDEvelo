Installation Guide
==================

This guide will walk you through the process of installing SDEvelo and its dependencies.

System Requirements
-------------------

Before installing SDEvelo, ensure your system meets the following requirements:

- **Operating Systems**: Linux (Ubuntu, CentOS), macOS, Windows 10
- **Python Version**: Python 3.8 and above
- **RAM**: Minimum 8GB, 16GB or more recommended for larger datasets
- **Storage**: At least 1GB of free disk space

Dependencies
------------

SDEvelo requires the following Python packages:

- anndata==0.10.7
- matplotlib==3.7.1
- numpy==1.23.5
- scipy==1.8.1
- scvelo==0.2.5
- seaborn==0.11.2
- torch==1.13.1+cu117

These will be automatically installed when you install SDEvelo.

Installation Steps
------------------

1. Set Up Python Environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We recommend using a virtual environment to avoid conflicts with other Python packages. You can create one using `venv` or `conda`:

Using venv:

.. code-block:: bash

   python3 -m venv sdevelo_env
   source sdevelo_env/bin/activate  # On Windows, use `sdevelo_env\Scripts\activate`

Using conda:

.. code-block:: bash

   conda create -n sdevelo_env python=3.8
   conda activate sdevelo_env

2. Install SDEvelo
^^^^^^^^^^^^^^^^^^

Once your environment is set up and activated, you can install SDEvelo using pip:

.. code-block:: bash

   pip install sdevelo

This command will install SDEvelo and all its dependencies.

3. Verify Installation
^^^^^^^^^^^^^^^^^^^^^^

To verify that SDEvelo has been installed correctly, you can run:

.. code-block:: bash

   python -c "import sdevelo; print(sdevelo.__version__)"

This should print the version number of SDEvelo without any errors.

Installing from Source
----------------------

If you want to install the latest development version of SDEvelo, you can install it directly from the GitHub repository:

.. code-block:: bash

   pip install git+https://github.com/Liao-Xu/SDEvelo.git

GPU Support
-----------

SDEvelo can leverage GPU acceleration for faster computations. If you have a CUDA-capable GPU, ensure you have the appropriate CUDA toolkit installed. The PyTorch version installed with SDEvelo (1.13.1+cu117) is compatible with CUDA 11.7.

To check if PyTorch can access your GPU, run:

.. code-block:: python

   import torch
   print(torch.cuda.is_available())

This should return `True` if your GPU is properly set up.

Troubleshooting
---------------

If you encounter any issues during installation:

1. Ensure you're using a supported Python version (3.8+).
2. Check that you have the latest version of pip: `pip install --upgrade pip`
3. If you're having issues with PyTorch, you may need to install it separately following the instructions on the `PyTorch website <https://pytorch.org/get-started/locally/>`_.
4. Make sure you have sufficient permissions to install packages on your system.

For any persistent issues, please refer to our GitHub issues page or contact our support team.

