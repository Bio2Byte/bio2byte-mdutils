#!/usr/bin/env python
#
# Author:   David Bickel
# Date:     14/12/2022

import argparse
import sys
from bio2byte.mdutils.pca import GmxDihedralPCA

def parse_arguments(arguments):
    
    parser = argparse.ArgumentParser("dihedral-pca", description="", add_help=False)

    inopt = parser.add_argument_group("Options to specify input files")
    inopt.add_argument("-s", "--struct", metavar="<.tpr/.gro/...>",
                       help="Structure file")
    inopt.add_argument("-f", "--traj", metavar="<.xtc/.trr/...>",
                       help="Trajectory file")

    outopt = parser.add_argument_group("Options to specify output files")
    outopt.add_argument("-o", "--eigenval", metavar="<.xvg>", default="eigenval.xvg",
                        help="Eigenvalues (default: eigenval.xvg)")
    outopt.add_argument("-v", "--eigenvec", metavar="<.trr>", default="eigenvec.trr",
                        help="Eigenvectors as full precision trajectory (default: eigenvec.trr)")
    outopt.add_argument("-j", "--proj", metavar="<.xvg>", default="proj.xvg",
                        help="Projection of trajectory along eigenvectors (default: proj.xvg)")
    outopt.add_argument("-n", "--index", metavar="<.ndx>", default="phi_psi.ndx",
                        help="Index file listing the dihedrals for the PCA in group 0 (default: phi_psi.ndx)")
    outopt.add_argument("--dihedral-trr", metavar="<.trr>", required=False,
                        help="Full precision trajectory to which the dihedrals are written. (Opt.)")
    outopt.add_argument("--dihedral-gro", metavar="<.gro>", required=False,
                        help="Structure file for the dihedral trajectory. (Opt.)")
    outopt.add_argument("--angdist", metavar="<.xvg>", required=False,
                        help="Data file to write the angle distribution to. (Opt.)")

    xopt = parser.add_argument_group("Other options")
    xopt.add_argument("--readv", action="store_true",
                      help="Do not calculate eigenvectors, but read them from provided file.")
    xopt.add_argument("--readn", action="store_true",
                      help="Do not generate dihedral indices, but read them  from provided file - Group 0")
    xopt.add_argument("--first", metavar="<int>", default="1",
                      help="First eigenvector for analysis (default: 1)")
    xopt.add_argument("--last", metavar="<int>", default="-1",
                      help="Last eigenvector for analysis (default: -1)")
    xopt.add_argument("-h", "--help", action="help")

    return parser.parse_args(arguments)


def main():
    args = parse_arguments(sys.argv[1:])
    pca = GmxDihedralPCA(args)
    pca.run()
    return 0

if __name__ == "__main__":
    sys.exit(main())
