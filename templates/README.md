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

### Template for ff99SB-type force fields

This template is used for force fields based on ff99SB, namely:
* ff99SB-ILDN
* a99SB-*disp*
* DES-Amber

### Template for CHARMM36m

This is a template for running simulations using the CHARM36m force field.
The following version is provided alongside.
* CHARMM36m (February 4, 2021) 


## Force fields

### AMBER ff99SB-ILDN  for GROMACS

Amber99SB-ILDN is an improvement to Amber99SB by Kresten Lindorff-Larsen and 
coworkers at D.E. Shaw Research, inc. [1,2]

The primary publication (please cite it) is:
Lindorff-Larsen, K., Piana, S., Palmo, K., Maragakis, P., Klepeis, J.L., 
Dror, R.O., Shaw, D.E., Improved side-chain torsion potentials for the Amber 
ff99SB protein force field, Proteins 78:1950–1958 (2010)

The present Gromacs implementation was taken from GROMCAS 2021-2. 

#### Added features

The following features have been added to the native force field:
* OPC water model [3]
* TIP4P/2005 water model [4]
* TIP4P/2005f water model [5]
* TIP4P-D water model [6]
* Parameters for phosphorylated residues [7,8]

### AMBER a99SB-*disp* for GROMACS

The a99SB-*disp* force field was developed to represent both folded and disordered 
protein states. It is based on ff99SB-ILDN with TIP4P-D water model. The vdW 
terms between water and amino acids were adjusted to improve the representation
of disordered protein states.[9]

The original version of the force field is available 
[here](https://github.com/paulrobustelli/Force-Fields).

#### Added features

The following features have been added to the native force field:
* Parameters for phosphorylated residues [7,8]

### DES-Amber for GROMACS

The a99SB-*disp* force field [9] revealed substential weaknesses modelling 
protein-protein interactions. As a result DES-Amber was developped, allowing for
more accurate simulations of protein−protein complexes, while still providing a 
state-of-the-art description of both ordered and disordered single-chain 
proteins.[10]

The original version of the force field is available 
[here](https://github.com/paulrobustelli/Force-Fields).

#### Added features

*None.*

### CHARMM36m for GROMACS

CHARMM36 force field in GROMACS format, including CGenFF version 4.4 and the 
CHARMM36m protein force field revision. Updated July 2020. Users can choose 
between C36 and C36m for protein simulations via the GROMACS "define" 
mechanism: `define = -DUSE_OLD_C36`.

The C36m parameter set is recommended for all protein simulations, but the 
ability to toggle between old and new parameter sets may be useful in the 
case of force field comparisons.[11]

The current CHARMM36 port for GROMACS is dated **February 4, 2021** and 
corresponds to the July 2020 toppar update.[12] Changes since March 2019 include
an update to CGenFF version 4.4, consolidation of parameters, new modified 
nucleic acid model compounds, and new NBFIXes for ions. The version posted in 
July 2020 omitted several NBFIXes related to Ca2+, which have been incorporated 
in the February 2021 version.

The original version of this force field is available 
[here](http://mackerell.umaryland.edu/charmm_ff.shtml#gromacs).

#### Added features

*None.*


## References

1.  Lindorff-Larsen, K. et al. (2010). Proteins 78, 1950-1958.
2.  Hornak, V. et al. (2006). Proteins 65, 712-725.
3.  Izadi, S. et al. (2014). J. Phys. Chem. Lett. 5, 3863-3871.
4.  Abascal, J.L.F. and Vega, C. (2005). J. Chem. Phys. 123, 234505.
5.  González, M.A. and Abascal, J.L.F. (2011). J. Chem. Phys. 135, 224516.
6.  Piana, S. et al. (2015). J. Phys. Chem. B 119, 5113-5123.
7.  Homeyer, N. et al. (2006). J. Mol. Model. 12, 281-289.
8.  Steinbrecher, T. et al. (2012). J. Chem. Theory Comput. 8, 4405-4412.
9.  Robustelli, P. et al. (2018). Proc. Natl. Acad. Sci. 115 (21), E4758-E4766.
10. Piana, S. et al. (2020). J. Chem. Theory Comput. 16 2494-2507.
11. Huang, J. et al. (2017) Nat. Methods, 14(1), 71–73. 
12. CHARMM port writted by: E. Prabhu Raman, Justin A. Lemkul, Robert Best, and Alexander D. MacKerell, Jr. ([source](http://mackerell.umaryland.edu/charmm_ff.shtml#gromacs), accessed 12/07/2021)
