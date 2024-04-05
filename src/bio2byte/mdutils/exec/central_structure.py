#!/usr/bin/env python
#
# Author: David Bickel 
# Date: 2022
"""
This module contains functionalities to approximate a centroid structure of a 
structural ensemble.

Rather than performing an all-vs-all RMSD fit (which is very time-consuming for
large structural ensembles), it applies a iterative approach:

    1. Pick a random initial structure from the ensemble as a reference structure
    2. Align the ensemble to the reference structure
    3. Calculate the average coordinates of the aligned structures
    4. Pick the structure with the lowest RMSD towards the average coordinates as new reference structure
    5. Repeat steps 2-4 until convergence is reached or the maximum number of iterations is exceeded

This methodology is likely to converge to a local centroid, rather than the absolute 
centroid. Therefore, the above steps are repeated multiple times to increase the
chance of finding a good approximation of the overall centroid.

help: ./central_structure.py -h
"""

import argparse
import os
import random
import subprocess
from subprocess import PIPE,STDOUT
import sys
import tempfile
import textwrap

import numpy as np

from bio2byte.mdutils.filetypes import GromacsXTC

os.environ["GMX_MAXBACKUP"] = "-1"


class gmxTrajectory:
    """ A class to store some general trajectory properties """

    def __init__(self, topology, trajectory):
        self.trajfile =   trajectory
        self.topfile  =   topology

        self.DEBUG    =   False
        self.LOG      =   []
        self._check_trajectory()

    def dumpFrame(self, dumpStructure, dumptime=-1.):
        if dumptime < 0:
            dumptime = random.uniform(self.t_start, self.t_stop)
        skip_to = dumptime - .6*self.t_step
        cmd = ["gmx", "trjconv",
               "-s", self.topfile, "-f", self.trajfile,
               "-b", "{:.1f}".format(skip_to),
               "-dump", "{:.1f}".format(dumptime),
               "-o", dumpStructure]
        g_trjconv = subprocess.run(cmd, input=b"System\n", capture_output=self.DEBUG)
        if self.DEBUG:
            self.LOG.append(g_trjconv)

    def avg_structure(self, avgStructure, fit_group="Backbone", covar_group="System", referenceStructure=None):
        if referenceStructure is None:
            referenceStructure = self.topfile
        tempXvg =  tempfile.mkstemp(suffix=".xvg", dir='.')[1]
        tempTrr =  tempfile.mkstemp(suffix=".trr", dir='.')[1]
        tempLog =  tempfile.mkstemp(suffix=".log", dir='.')[1]

        cmd = ["gmx", "covar",
               "-s", referenceStructure, "-f", self.trajfile,
               "-o", tempXvg, "-v", tempTrr, "-l", tempLog,
               "-av", avgStructure]
        cmd_input = (fit_group+'\n'+covar_group+'\n').encode()
        g_covar = subprocess.run(cmd, input=cmd_input, capture_output=self.DEBUG)
        if self.DEBUG:
            self.LOG.append(g_covar)
        os.remove(tempXvg)
        os.remove(tempTrr)
        os.remove(tempLog)

    def rmsd(self, referenceStructure=None, fit_group="C-alpha", rms_group="C-alpha"):
        if referenceStructure is None:
            referenceStructure = self.topfile
        tempXvg =  tempfile.mkstemp(suffix=".xvg", dir='.')[1]

        cmd = ["gmx", "rms",
               "-s", referenceStructure, "-f", self.trajfile,
               "-o", tempXvg, "-xvg", "none"]
        cmd_input = (fit_group+'\n'+rms_group+'\n').encode()
        g_rms   = subprocess.run(cmd, input=cmd_input,
                                 capture_output=self.DEBUG)
        if self.DEBUG:
            self.LOG.append(g_rms)
        rmsdarr = np.genfromtxt(tempXvg, dtype=np.float32)
        os.remove(tempXvg)
        return rmsdarr

    def _check_trajectory(self):
        xtc = GromacsXTC(self.trajfile)
        self.t_start = xtc.time0
        self.t_stop  = xtc.time1
        self.frames  = xtc.nframes
        self.t_step  = (self.t_stop - self.t_start) / (xtc.nframes - 1)


class centroidFinder:
    """ Class that helps finding centroids iteratively """

    def __init__(self, topology_file, trajectory_file):
        self.topology_file   =  topology_file
        self.trajectory_file =  trajectory_file
        self.gmxTraj         =  gmxTrajectory(topology_file,
                                              trajectory_file)
        # Centroid suggestions
        self.centroids = dict()

    def __call__(self, outfile, repeats=10, maxiter=10):
        for n in range(repeats):
            self.find_centroid(maxiter=maxiter)

        for idx,(pdb, disp) in self.centroids.items():
            try:
                assert mindisp <= disp
                os.remove(pdb)
            except UnboundLocalError:
                minidx  = idx
                mindisp = disp
                minpdb  = pdb
            except AssertionError:
                os.remove(minpdb)
                minidx  = idx
                mindisp = disp
                minpdb  = pdb
            #except FileNotFoundError:
                #pass
        print(self.centroids)
        print("rename: %s -> %s" % (minpdb, outfile))
        os.rename(minpdb, outfile)
        sys.stdout.write("Centroid structure: {0}\n".format(outfile))
        sys.stdout.write("  Frame: {:d}\n".format(minidx+1))
        sys.stdout.write("  D^2:   {:.6f}\n".format(mindisp))

    def find_centroid(self, maxiter=5,):
        referencePDB = tempfile.mkstemp(suffix=".pdb", dir='.')[1]
        averagePDB   = tempfile.mkstemp(suffix=".pdb", dir='.')[1]
        k = -1
        while k in self.centroids:
            k -= 1
        self.centroids[k] = [referencePDB, np.inf]

        self.gmxTraj.dumpFrame(referencePDB)

        for i in range(maxiter):
            #sys.stderr.write("Iteration {0:d}".format(i))

            # Calculate average structure
            self.gmxTraj.avg_structure(averagePDB, fit_group="Backbone", covar_group="System",
                                       referenceStructure=referencePDB)
            #sys.stderr.write('.')

            # RMSD-fit to average structure
            rms_arr = self.gmxTraj.rmsd(averagePDB, fit_group="Backbone", rms_group="Backbone")
            min_idx = np.argmin(rms_arr[:,1])
            min_t   = rms_arr[min_idx,0]
            if min_idx in self.centroids:
                #sys.stderr.write(".\n")
                break
            #sys.stderr.write('.')

            # Get closest to average structure
            referencePDB = tempfile.mkstemp(suffix=".pdb", dir='.')[1]
            self.gmxTraj.dumpFrame(referencePDB, dumptime=min_t)
            sys.stderr.write('.')

            # Evaluate new dispersion
            rms_arr = self.gmxTraj.rmsd(referencePDB, fit_group="Backbone", rms_group="Backbone")
            rmsd    = rms_arr[:,1]
            dispersion = np.sum(rmsd**2) / (rmsd.size - 1)
            self.centroids[min_idx] = [referencePDB, dispersion]
            sys.stderr.write('.\n')
        else:
            pass  #sys.stderr.write("[ WARNING ] Iterative averaging did not converge.\n")

        referencePDB,dispersion = self.centroids[min_idx]
        #sys.stderr.write("[ {0:d} ]\n".format(min_idx))
        #sys.stderr.write("  File:\t{0}\n".format(referencePDB))
        #sys.stderr.write("  D^2:\t{0:.6f}\n".format(dispersion))
        os.remove(averagePDB)
        return min_idx


def parse_cmdline(cmdline_arguments):
    """ Parse command line parameters """
    AParser =   argparse.ArgumentParser(description=textwrap.dedent("""\
        This script contains functionalities to approximate a centroid structure of a 
        structural ensemble.

        Rather than performing an all-vs-all RMSD fit (which is very time-consuming for
        large structural ensembles), it applies a iterative approach:

            1. Pick a random initial structure from the ensemble as a reference structure
            2. Align the ensemble to the reference structure
            3. Calculate the average coordinates of the aligned structures
            4. Pick the structure with the lowest RMSD towards the average coordinates as new reference structure
            5. Repeat steps 2-4 until convergence is reached or the maximum number of iterations is exceeded
        """))

    arg_io  =   AParser.add_argument_group('Input/output options')
    arg_io.add_argument("-s", "--struct", metavar="<GRO>", required=True,
        help=("Structural information for the coordinates in the trajectory: [.tpr, .gro, .pdb]"))
    arg_io.add_argument("-f", "--traj", metavar="<XTC>", required=True,
        help="Trajectory file: [.trr, .xtc]")
    arg_io.add_argument("-o", "--outfile", metavar="<PDB>", default="centroid.pdb",
        help="Output file with centroid structure")

    arg_opt  =   AParser.add_argument_group('Options')
    arg_opt.add_argument("-R", "--repeats", default=10, type=int, metavar="<INT>",
        help="Sets how many independent repeats are done to identify the global minimum. [default: 10]")
    arg_opt.add_argument("--maxiter", default=5, type=int, metavar="<INT>",
        help="Defines the maximal number of local minimization steps before a local minimum is reported. [default: 5]")

    args  =  AParser.parse_args(cmdline_arguments)
    if args.repeats < 1:
        raise ValueError("[ERROR] Too few repeats: %d; must be at least 1" % args.repeats)
    if args.maxiter < 2:
        raise ValueError("[ERROR] Too small maxiter: %d; must be at least 2 (recommended 4-5)." % args.maxiter)
    return args


def main():
    # Parse arguments
    args = parse_cmdline(sys.argv[1:])
    # Find centroid
    CFind = centroidFinder(args.struct, args.traj)
    CFind(args.outfile, repeats=args.repeats, maxiter=args.maxiter)
    return 0


if __name__ == "__main__":
    sys.exit(main())
