# MD Template Structure #

This repository contains a MD template structure, that can be easily modified 
to run arbitrary simulations. It contains scripts to facilitate the system
preparation and submission as well as a folder structure to structure the 
obtained data.

## Supported force fields

* ff99SB-ILDN [1]
* a99SB-disp [2]
* DES-Amber [3]

## Folder structure ##

```
    .
    │
    ├── 1_prepi
    │   └── preparation.sh
    │
    ├── 2_equil
    │   ├── 01_min_wat.mdp
    │   ├── 02_min_slt.mdp
    │   ├── 03_min_all.mdp
    │   ├── 04_nvt_heat.mdp
    │   ├── 05_npt_equi.mdp
    │   └── equilibration.slrm
    │
    └── 3_prod
        ├── npt_pme.mdp
        └── production.slrm
```

### ./1_prepi

Preparation of all input structures.
The script `preparation.sh` takes a PDB-file and processes it interactively:

1.  Prepare PDB: `pdb2gmx`
2.  Create pbc box: `editconf`
3.  Solvate: `solvate`
4.  Add ions: `genion`
5.  Generate a snapshot of the prepared system

Final files generated:

*   Initial configuration: **COORD.gro**
*   Topology: **TOPOL.top**
*   Snapshot (pdb) of the prepared system: **COORD.gro.pdb**

### ./2_equil

Here, equilibration MD is done. Any number of mdp-files can be added.
The `equilibration.sh` will process these files in an alpha-numerical order.
Comment lines at the top of each MDP-file should provide a short description of
the equilibration step.

Standard files:

1.  `01_min_wat.mdp`:   Minimization of solvent
2.  `02_min_slt.mdp`:   Minimization of solute
3.  `03_min_all.mdp`:   Minimization of full system
4.  `04_nvt_heat.mdp`:  Temperature relaxation by short NVT MD
5.  `05_npt_equi.mdp`:  Density relaxation by short NPT MD

### ./3_prod

Here, the production MD is executed using NPT ensemble. The file `npt_pme.mdp` 
provides the simulation parameters.

The SLURM script `production.slrm` cn be used to run simulations on  a HPC 
cluster (e.g., Hydra). By changing the `N_ITERATIONS` variable, the desired 
simulation time can be adjusted.

## References

1. Lindorff-Larsen, K. et al. (2010). Proteins 78, 1950-1958.
2. Robustelli, P. et al. (2018). PNAS 115 (21), E4758-E4766.
3. Piana, S. et al. (2020). JCTC, 16 2494-2507.
