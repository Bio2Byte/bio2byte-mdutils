#!/usr/bin/env python 

import argparse
import os
import sys
import numpy as np
from bio2byte.mdutils.utils.ndhist import NDHistogram


def parse_arguments(arguments):
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "-i", "--infile",)
    parser.add_argument("-o", "--outfile", default="")
    parser.add_argument("--frame-window", type=int, default=500)
    args = parser.parse_args(arguments)

    if args.outfile == "":
        if os.path.dirname(args.infile) != "":
            args.outfile = f"{os.path.dirname(args.infile)}/replica_correff8.dat"
        else:
            args.outfile = "replica_correff8.dat"
    
    return args


def correlation_coeff(hists, n_samples):
    arr = hists / n_samples
    corr = []
    for i, arri in enumerate(arr):
        mi = (arri > 0)
        for arrj in arr[i+1:]:
            mj = (arrj > 0)
            ai, aj = arri[mi | mj], arrj[mi | mj]
            sp = np.sum( (ai - np.mean(ai)) * (aj - np.mean(aj)) )
            sqi, sqj = np.sum((ai - np.mean(ai))**2), np.sum((aj - np.mean(aj))**2)
            corr.append( sp / np.sqrt(sqi * sqj) )
    return corr


def main():

    # Parse arguments
    args = parse_arguments(sys.argv[1:])

    # Read projections
    proj = NDHistogram.from_projectionfile(args.infile, n_replicas=4)
    proj.data = proj.data[:,:,:8]
    proj.make_dbins(nbins=5)

    with open(args.outfile, "w") as fhandle:
        fhandle.write(
            "#NFrames\t" + \
            "\t".join([f"r{i+1}_r{j+1}" for i in range(proj.n_replicas) for j in range(i+1, proj.n_replicas)]) + "\n"
        )
        for n_frames, hists in proj.iter_histograms(first=0, last=args.frame_window, inc_last=args.frame_window):
            corr = correlation_coeff(hists, n_frames)
            fhandle.write("{0}\t{1}\n".format(
                n_frames, '\t'.join([f'{v:.4f}' for v in corr])
            ))

if __name__ == "__main__":
    sys.exit(main())
