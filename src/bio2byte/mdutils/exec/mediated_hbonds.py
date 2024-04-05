#!/usr/bin/env python
import argparse
import sys
import textwrap as tw

import MDAnalysis as mda
from bio2byte.mdutils.hbonds import MediatedHBondAnalyzer

def parse_cmdline(cmdline_arguments):
    """ Parse command line parameters """

    aparser = argparse.ArgumentParser(description=tw.dedent("""\
        Script to analyze (solvent-mediated) hydrogen bonds in a md simulations
        """), add_help=False)
    
    # Input/output
    ioargs = aparser.add_argument_group('Input/output options')
    ioargs.add_argument("--structure", "-s", 
        help=("Structure containing topology information for the trajectory."
              "[tpr, pdb, gro, top]"))
    ioargs.add_argument("--trajectory", "-f", nargs="+",
        help=("MD trajectory containing the coordinates to analyze. "
              "[xtc, trr, nc, mdcrd]"))
    ioargs.add_argument("--outfile", "-o", default="hbonds.csv",
        help=("(Optional) File to write the output to. "
              "Default: hbonds.csv"))
    
    selargs = aparser.add_argument_group('Selection options')
    selargs.add_argument("--sele0", default="all", required=False,
        help=("(Optional) If defined, only h-bonds within sele0 are calculated. "
              "If in addition sele1 is defined, only h-bonds between sele0 and "
              "sele1 are calculated."), metavar="SELECTION")
    selargs.add_argument("--sele1", default=None, required=False,
        help=("(Optional) If sele1 is defined, only hydrogen bonds between "
              "sele0 and sele1 are calculated."), metavar="SELECTION")
    selargs.add_argument("--solvent", default="", required=False,
        help=("(Optional) If defined h-bonds between sele0 and sele1 can be "
              "mediated by a single solvent molecule."), metavar="SELECTION")
    selargs.add_argument("--update-sele", action="store_true",
        help="(Optional) Update the selection dynamically for each frame.")

    optargs = aparser.add_argument_group('Options')
    optargs.add_argument("-d", "--distance-cutoff", type=float, default=3.0,
        help=("(Optional) Defines the maximal donor-acceptor distance in "
              "Angstrom to be considered a hydrogen bonds . Default: 3.0"))
    optargs.add_argument("-a", "--angle-cutoff", type=float, default=150.0,
        help=("(Optional) Defines minimum angle donor-H-acceptor angle in"
              "degrees to be considered a hydrogen bond. Default: 150.0"))
    optargs.add_argument("--append", action="store_true", default=False,
        help=("(Optional) If this argument is given, the results will be "
              "appended to an existing output file, rather than overwriting it"))
    aparser.add_argument("--help", "-h", action="help", 
        help="Shows this help message.")
    
    args = aparser.parse_args(cmdline_arguments)
    return args

def main():
    # Parse command line arguments
    args = parse_cmdline(sys.argv[1:])

    # Initialize universe
    MD = mda.Universe(args.structure, args.trajectory)

    # Run hydrogen bond analysis
    analyzer = MediatedHBondAnalyzer(
        MD, args.sele0, args.sele1, args.solvent, 
        args.distance_cutoff, args.angle_cutoff, args.update_sele)
    analyzer.run()

    # Write results
    if args.append:
        handler = open(args.outfile, "a")
    else:
        handler = open(args.outfile, "w")
        handler.write("#frame,resn0,resid0,resn1,resid1,type\n")
    for hb in analyzer.results:
        handler.write(hb.to_csv() + "\n")
    handler.close()

if __name__ == "__main__":
    sys.exit(main())