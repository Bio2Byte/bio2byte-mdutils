#!/usr/bin/env python

import sys
import argparse

from bio2byte.mdutils.dihedrals import DihedralAnalyzer

def parse_cmdline(cmdline_arguments):
    """ Parse command line parameters """
    AParser =   argparse.ArgumentParser(description="")

    arg_io  =   AParser.add_argument_group('Input/output options')
    arg_io.add_argument("-s", "--struct",
        metavar="<GRO>", required=True,
        help=("Structural information for the coordinates "
              "in the trajectory: [.tpr, .gro, .pdb]"))
    arg_io.add_argument("-f", "--traj",
        nargs='+', metavar="<XTC>", required=True,
        help="Trajectory file(s): [.trr, .xtc]")
    arg_io.add_argument("-o", "--outfile",
        metavar="<CSV>", default="dihedrals.csv",
        help="Output datafile with dihedral data.")

    arg_opt  =   AParser.add_argument_group('Options')
    arg_opt.add_argument("-D", "--dihedrals",
                         default="phi,psi",
                         help=("Comma-separated list of dihedrals to calculate:"
                               " phi, psi, chi1. [default: \"phi,psi\"]"))
    arg_opt.add_argument("--mask", default="all")
    arg_opt.add_argument("--unwrap", action="store_true")
    arg_opt.add_argument("--precision", default=2, type=int,
                         help=("Determines the output precision for the "
                               "dihedral angles. [default: 2]"))
    # --otype [csv,dat,pkl]
    # --delim
    
    args  =  AParser.parse_args(cmdline_arguments)
    return args


def main():
    # Parse command line arguments
    args = parse_cmdline(sys.argv[1:])

    # Initialize DihedralAnalyzer
    DA = DihedralAnalyzer(args)

    # Calculate dihedrals
    dihedrals = args.dihedrals.lower().split(",")
    dihedral_list = [d.strip() for d in dihedrals]

    dihedralDF = DA.calculate_dihedrals(dihedral_list)

    dihedralDF.to_csv(
        args.outfile, sep=',', na_rep="nan", 
        float_format="%.{0:d}f".format(args.precision),
        index_label="#Frame")

    return 0



if __name__ == "__main__":
    sys.exit(main())
