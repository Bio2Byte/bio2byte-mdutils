#!/usr/bin/env python
import argparse
import os, sys
from tempfile import NamedTemporaryFile
import numpy as np
import pandas as pd

from constava.__main__ import parse_parameters, run_analyze


def convert_dataframe(input_dataframe):
    outdf = pd.DataFrame({
        "#Frame":   [],
        "ResIndex": [],
        "ResName":  [],
        "Phi[rad]": [],
        "Psi[rad]": []
    })
    for phicol in input_dataframe.columns:
        if not phicol.endswith(":phi"):
            continue
        resn, resi, _ = phicol.split(":")
        psicol = f"{resn}:{resi}:psi"
        resdf = input_dataframe[["#Frame", phicol, psicol]]
        resdf = resdf.rename(columns={phicol: "Phi[rad]", psicol: "Psi[rad]"})
        resdf["ResIndex"] = resi
        resdf["ResName"] = resn
        outdf = pd.concat([outdf, resdf], axis=0)
    outdf = outdf.astype({"#Frame": int})
    outdf["Phi[rad]"] = np.radians(outdf["Phi[rad]"])
    outdf["Psi[rad]"] = np.radians(outdf["Psi[rad]"])
    return outdf


def main():
    # Parse command line arguments
    args = parse_parameters(sys.argv[1:])

    # Read and convert the dihedral dataframe to constava format
    df = pd.read_csv(args.input_file)
    df = convert_dataframe(df)
    with NamedTemporaryFile(dir=os.getcwd(), suffix=".csv", delete=True) as tmpcsv:
        df.to_csv(tmpcsv, index=False)
    
    # Parse command line parameters
    args = parse_parameters(sys.argv[1:])
    assert args.subcommand == "analyze", "ERROR: Can only run 'analyze' in batch"

    # Convert the input csv to constava format
    df = pd.read_csv(args.input)
    df = convert_dataframe(df)
    with NamedTemporaryFile(dir=os.getcwd(), suffix=".csv", delete=True) as tmpcsv:
        df.to_csv(tmpcsv, index=False)
        args.input = tmpcsv.name
        run_analyze(args)

if __name__ == "__main__":
    sys.exit(main())