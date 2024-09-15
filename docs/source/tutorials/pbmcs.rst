# Peripheral Blood Mononuclear Cells (PBMCs) Tutorial

This tutorial demonstrates the application of SDEvelo to Peripheral Blood Mononuclear Cells (PBMCs) data. PBMCs include a diverse array of mature immune cells, making them an excellent model for studying immune system dynamics and cellular heterogeneity. Importantly, this tutorial also serves as a negative control for transcriptional dynamics in mature-state cell populations.

## Background

PBMCs in their mature state are expected to lack dynamic information, as the mRNA levels of these cells have already equilibrated. This makes them an ideal dataset for demonstrating SDEvelo's ability to accurately detect the absence of directional patterns in steady-state populations, in contrast to other methods that may erroneously infer strong directional patterns.

## Dataset

The analysis is performed on a PBMC dataset consisting of 65,877 cells and 33,939 genes, generated using the 10x platform. After quality control (QC), RNA velocity analysis was applied to the remaining 601 genes.

## Contents

In this tutorial, we will:

1. Load and preprocess the PBMC single-cell RNA sequencing data
2. Apply SDEvelo to model the transcriptional dynamics in different immune cell populations


## Key Findings

- SDEvelo detects random directions between most cell types, consistent with the expected behavior of mature-state populations.
- The inferred latent time from SDEvelo shows most cells close to the end of their differentiation, consistent with the terminal cell types seen among PBMCs.

## Notebook

.. toctree::
   :maxdepth: 1

   demo_pbmc
