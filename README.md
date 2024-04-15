# Bio2byte :: Molecular dynamics utilities

*by David Bickel and Wim Vranken*

This repository contains an assembly of templates, scripts, and Jupyter 
notebooks for running and analyzing large-scale MD simulations.

## Contents

```
├── notebooks/                              # Jupyter notebook templates
│   ├── bmrb_comparison.ipynb
│   ├── build_peptides.ipynb
│   ├── constava_adjacent.ipynb
│   ├── constava.ipynb
│   ├── dihedral_timeseries.ipynb
│   ├── dssp_summary.ipynb
│   ├── hbond_plots.ipynb
│   ├── pca_comparison.ipynb
│   ├── pca_replica_correff.ipynb
│   ├── ramachandran.ipynb
│   ├── rama_effect_propagation.ipynb
│   ├── rama_rmsd.ipynb
│   └── README.md
├── slurm/                                  # SLURM script templates
│   ├── conformational_variability.slrm
│   ├── dihedral_pca.slrm
│   ├── dihedrals.slrm
│   ├── mediated_hbonds.slrm
│   ├── prepare_mdsystems.slrm
│   ├── process_trajectories.slrm
│   ├── replica_correlation.slrm
│   └── README.md
├── templates/                              # MD simulation templates
│   ├── forcefields/
│   ├── template_ff99sb/
│   └── README.md
├── src/                                    # Python code for analyses
│   └── bio2byte/mdutils
└── README.md (this file)
```

## Building and installing the repository

If you just want to use the templates or copy some code from the notebooks,
there is no need to install the package. Just go through the directories and 
copy what you need.

If you though, want to use the prewritten scripts, SLURM scripts and notebooks
as they are, you will need to build this package and install it. Here, is how:

```sh
# Create a virtual environment (optional but recommended)
python -m venv mdutils
source mdutils/bin/activate

# Get a copy of the repository on your local machine
git clone https://github.com/Bio2Byte/bio2byte-mdutils.git

# Enter the repository
cd bio2byte-mdutils

# Build the repository
make build

# Install the repository
make install
```

Installing the repository will give you access to a couple of commands to help
you with your MD analyses (also see Usage):

| Command                       | Description |
|-------------------------------|-------------|
| `mdutils-process-trajectory`  | Processes MD trajectories (center protein, PBC image, remove solvent) |
| `mdutils-hbonds-wrapper`      | Calculates solute hydrogen bonds as well as solvent mediated hydrogen bonds |
| `mdutils-hbonds`              | Calculates hydrogen bonds in a MD trajectory |
| `mdutils-centroid`            | Estimates the centroid structure of a MD ensemble |
| `mdutils-dihedral`            | Calculates phi/psi dihedral angels in a MD ensemble |
| `mdutils-dihedral-pca`        | Performes a PCA on the phi/psi backbone dihedrals |
| `mdutils-replica-correlation` | Calculates correlation coefficients between replicas along N principal components |
| `mdutils-ionic-contacts`      | Calculates "ionic contacts" in the MD ensemble |
| `mdutils-protein-contacts`    | Calculates protein contacts in the MD ensemble |
| `mdutils-constava-wrapper`    | Infers residue-wise conformational states from the MD ensemble |

## Citation

The scripts in this repository were generated over the course of a 
large-scale simulation project. When citing this repository, please
cite the associated [publication](https://www.biorxiv.org/content/10.1101/2024.02.15.580491v1):

> Bickel, D. & Vranken, W. F. 
> Effects of phosphorylation on protein backbone dynamics and conformational preferences. 
> bioRxiv (2024) doi:10.1101/2024.02.15.580491.

## License

This repository is licensed under MIT license.

## Contact

[![ORCID](https://orcid.org/sites/default/files/images/orcid_16x16.png)](https://orcid.org/0000-0003-0332-8338) - David Bickel - [david.bickel@vub.be](mailto:david.bickel@vub.be)

[![ORCID](https://orcid.org/sites/default/files/images/orcid_16x16.png)](https://orcid.org/0000-0001-7470-4324) - Wim Vranken -  [wim.vranken@vub.be](mailto:wim.vranken@vub.be)
