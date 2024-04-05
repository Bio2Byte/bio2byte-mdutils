# SLURM submission scripts notebooks

This folder contains notebooks that were used for performing analyses of the 
MD simulations and generating figures for publication.

## Overview

| File | Description |
|------|-------------|
| **Preparation** | |
| `prepare_mdsystems.slrm` | Run the preparation and equilibration in the MD directory |
| **Analysis** | |
| `process_trajectories.slrm` | Process MD trajectories (center on protein, image, and remove solvent molecules) |
| `dihedrals.slrm` | Calculate dihedral angles along a MD trajectory |
| `conformational_variability.slrm` | Calculate the conformational state variablility of MD ensembles |
| `mediated_hbonds.slrm` | Calculate hydrogen bonds along MD trajectories |
| `dihedral_pca.slrm` | Perform a backbone dihedral PCA on the trajectory |
| `replica_correlation.slrm` | Calculate the correlation coefficients between replicas in principal component-space |
