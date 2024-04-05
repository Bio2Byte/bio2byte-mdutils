#!/usr/bin/env python
import sys
import argparse
import textwrap
from multiprocessing import Pool

import MDAnalysis as mda

from bio2byte.mdutils.contacts import ContactAnalyzer
from bio2byte.mdutils.filetypes import GromacsNDX


def parse_arguments(arguments):
    """ Parse command line arguments passed by th user """

    aparser = argparse.ArgumentParser(description=textwrap.dedent("""\
        Script to analyze contacts between two selections in a md simulation.
        """))

    # Input/output
    aparser.add_argument("--structure", "-s", 
        help=("Structure containing topology information for the trajectory."
              "[tpr, pdb, gro, top]"))
    aparser.add_argument("--trajectory", "-f", nargs="+",
        help=("MD trajectory containing the coordinates to analyze. "
              "[xtc, trr, nc, mdcrd]"))
    aparser.add_argument("--outfile", "-o", default="contacts.csv",
        help=("(Optional) File to write the output to. "
              "Default: contacts.csv"))
    aparser.add_argument("--ndx", "-n", default=None,
        help=("(Optional) Index file for selecting groups in a GROMACS-style "
              "manner. This overrules selections passed by --atom-selection1 "
              "and --atom-selection2. Default: none"))
    
    # Options
    aparser.add_argument("--atom-selection1", "-a", default="name CB",
        help=("(Optional) Atomic selection(s) between which the ionic contacts"
              "are calculated. Default: 'Protein'"))
    aparser.add_argument("--atom-selection2", "-b", default="name CB",
        help=("(Optional) Atomic selection(s) between which the ionic contacts"
              "are calculated."))
    aparser.add_argument("--distance-cutoff", "-c", type=float, default=6.,
        help=("(Optional) Maximum distance in Angstrom between two atoms, to still "
              "count as contacts. Default: 6."))
    aparser.add_argument("--byres", action="store_true",
        help=("(Optional) Only report the lowest distance contacts."))
    aparser.add_argument("--nprocs", "-j", type=int, default=1,
        help=("(Optional) Set to value .gt. 1 to allow multicore processing. "
              "Default: 1"))
    
    args = aparser.parse_args(arguments)
    return args


def main():
    # Parse commandline arguments
    args = parse_arguments(sys.argv[1:])

    # Load Universe
    MD = mda.Universe(args.structure, args.trajectory)

    # Get selections
    if args.ndx is None:
        sele0 = MD.select_atoms(args.atom_selection1)
        sele1 = MD.select_atoms(args.atom_selection2)
    else:
        sele0, sele1 = None, None
        ndx = GromacsNDX(args.ndx, index_correction=-1)
        print("\nSelect two groups, between which contacts will be calculated.\n")
        for i,lbl in enumerate(ndx._group_labels):
            print("{0:8d}:  {1:16}  ({2:d} atoms)".format(i, lbl, ndx[lbl].size))
        print()

        strinp = input("Select first group: ")
        if strinp.isdigit(): strinp = int(strinp)
        sele0 = MD.atoms[ndx[strinp]]
        print("Group with {0:d} atoms selected.\n".format(sele0.n_atoms))

        strinp = input("Select second group: ")
        if strinp.isdigit(): strinp = int(strinp)
        sele1 = MD.atoms[ndx[strinp]]
        print("Group with {0:d} atoms selected.\n".format(sele1.n_atoms))
 
    # Contact analysis
    Contacts = ContactAnalyzer(MD, sele0, sele1,
                distance_cutoff = args.distance_cutoff,
                by_residue = args.byres,
                nprocs = args.nprocs)
    Contacts.run()

    # Write resutls
    with open(args.outfile, "w") as f:
        f.write("#time,resid_0,atomid_0,resid_1,atomid_1,distance\n")
        for cnts in Contacts.result:
            f.write("{0:.1f},{1:d},{2:d},{3:d},{4:d},{5:.3f}\n".format(*cnts))

    return 0


if __name__ == "__main__":
    sys.exit(main())
