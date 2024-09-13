SDEvelo Documentation
=====================

.. image:: _static/sde_flow.png
   :alt: SDEvelo Flow Chart
   :align: center

SDEvelo: A Deep Generative Approach for Transcriptional Dynamics
----------------------------------------------------------------

SDEvelo represents a significant advancement in single-cell RNA sequencing (scRNA-seq) an sequece based spatial transcriptomics (ST) analysis:

-  SDEvelo employs multivariate stochastic differential equations (SDE) to model RNA velocity.

- **Key Features**:

  - Deep generative approach
  - Models dynamics of unspliced and spliced RNAs
  - Captures inherent uncertainty in transcriptional dynamics
  - Estimates cell-specific latent time across genes

SDEvelo offers a more accurate and comprehensive approach to understanding cell differentiation and state transitions in scRNA-seq studies.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   overview
   installation



Key Features
------------

- Multivariate stochastic modeling of transcriptional dynamics
- Cell-specific latent time estimation
- Applicable to both scRNA-seq and sequencing-based spatial transcriptomics data
- Computationally scalable
- Accurate detection of carcinogenesis
- Facilitates downstream analyses for biological discovery

System Requirements
-------------------

- **Operating Systems**: Linux (Ubuntu, CentOS), macOS, Windows 10
- **Python Version**: Python 3.8 and above
- **Hardware Requirements**: No non-standard hardware required
- **Installation Time**: Approximately 5 minutes

For more detailed information, please refer to the section.

Quick Start
-----------

.. code-block:: bash

   pip install sdevelo

For more detailed usage instructions, please see the  section.

Demo
----

We provide an interactive demo as a Jupyter Notebook. For details on running the demo and expected outputs, please check the section.



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`