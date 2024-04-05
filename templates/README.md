# MD templates

The repository contains templates and instructions for running simulations with
common state-of-the-art force fields. Some of these force fields have added 
features (e.g., phosphate residues, newer water models).

The `forcefields/` folder contains all force fields while the `template_*/`
directories contain templates to run simulations with those force fields.


## Templates

All template directories follow a common structure, e.g.:

```
template_ff99sb/
├── 1_prepi
│   ├── ff99sb-ildn-phosaa.ff -> ../../forcefields/ff99sb-ildn-phosaa.ff
│   └── preparation.sh
├── 2_equil
│   ├── 01_min_wat.mdp
│   ├── 02_min_slt.mdp
│   ├── 03_min_all.mdp
│   ├── 04_nvt_heat.mdp
│   ├── 05_npt_equi.mdp
│   └── equilibration.slrm
├── 3_prod
│   ├── npt_pme.mdp
│   └── production.slrm
└── README.md
```

In the `1_prepi/` directory, the input files for the simulation are prepared. 
This can be done by running `./preparation.sh <INPUT_STRUCTURE.pdb>` 
inside the directory

In the `2_equil/` directory an initial equilibration of the MD system 
is performed, consisting of a three-step energy minimization 
(`01_min_wat.mdp`, `02_min_slt.mdp`, `03_min_all.mdp`), a temperature
equilibration (`04_nvt_heat.mdp`), and a pressure equilibration 
(`05_npt_equi.mdp`). The SLURM submission script `equilibration.slrm`
will run all *.mdp files in order of their appearance. 

In the `3_prod/` finally contains the production simulation. The
`production.slrm` will run one iteration (by default 100 ns) of the
simulation, and resubmit itself until the maximum number of iterations
is reached.

## Force fields

## AMBER ff99SB-ILDN  for GROMACS

Amber99SB-ILDN is an improvement to Amber99SB by Kresten Lindorff-Larsen and 
coworkers at D.E. Shaw Research, inc. [1,2]

The primary publication (please cite it) is:
Lindorff-Larsen, K., Piana, S., Palmo, K., Maragakis, P., Klepeis, J.L., 
Dror, R.O., Shaw, D.E., Improved side-chain torsion potentials for the Amber 
ff99SB protein force field, Proteins 78:1950–1958 (2010)

The present Gromacs implementation was taken from GROMCAS 2021-2. 

### Added features

The following features have been added to the native force field:
* OPC water model [3]
* TIP4P/2005 water model [4]
* TIP4P/2005f water model [5]
* TIP4P-D water model [6]
* Parameters for phosphorylated residues [7,8]

### References

1.  Lindorff-Larsen, K. et al. (2010). Proteins 78, 1950-1958.
2.  Hornak, V. et al. (2006). Proteins 65, 712-725.
3.  Izadi, S. et al. (2014). J. Phys. Chem. Lett. 5, 3863-3871.
4.  Abascal, J.L.F. and Vega, C. (2005). J. Chem. Phys. 123, 234505.
5.  González, M.A. and Abascal, J.L.F. (2011). J. Chem. Phys. 135, 224516.
6.  Piana, S. et al. (2015). J. Phys. Chem. B 119, 5113-5123.
7.  Homeyer, N. et al. (2006). J. Mol. Model. 12, 281-289.
8.  Steinbrecher, T. et al. (2012). J. Chem. Theory Comput. 8, 4405-4412.
