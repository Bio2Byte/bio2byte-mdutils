# Jupyter notebooks

This folder contains notebooks that were used for performing analyses of the 
MD simulations and generating figures for publication.

## Overview

| File                          | Description |
|-------------------------------|---|
|  **MD Preparation**           | |
| `build_peptides.ipynb`        | This script was used to build the starting conformations of the peptides |
| **Analyses**                  | |
| `dihedral_timeseries.ipynb`   | Calculates circular variance onf phi/psi dihedrals |
| `ramachandran.ipynb`          | Visualization of dihedral angles as ramachandran plots |
| `ramachandran_rmsd.ipynb`     | Comparison between peptides, based on their ramachandran profiles |
| `rama_effect_propagation.ipynb` | Visualization of how residues affect sampling of adjacent residues |
| `dssp_summary.ipynb`          | Analyzes and visualizes the DSSP results fro the simulations |
| `constava.ipynb`              | Visualizes the conformational state propensities for all residues |
| `constava_adjacent.ipynb`     | Visualizes the effects of residues on the conformational state propensities of their neighbors |
| `hbond_plots.ipynb`           | Visualizes the hydrogen bond networks in the peptides |
| `pca_comparison.ipynb`        | Comparative nalysis and visualization of PCA results between peptides |
| `pca_replica_correff.ipynb`   | Visualization of  PCA results and time series of the replica correlation |
| `bmrb_comparison.ipynb`       | Comparison of the simulation results with NMR data from the [BMRB](https://bmrb.io/) |