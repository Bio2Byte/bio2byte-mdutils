#!/usr/bin/env python
import os, sys
import numpy as np
import pandas as pd


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

def main(infile, outfile = None):
    if outfile is None:
        fpath,fext = os.path.splitext(infile)
        outfile = f"{fpath}2{fext}"
    # Load dataframe
    input_dataframe = pd.read_csv(infile)
    # Convert
    output_dataframe = convert_dataframe(input_dataframe)
    # Write output dataframe
    output_dataframe.to_csv(outfile, float_format="%.5f", index=False)


if  __name__ == "__main__":
    with pd.option_context("mode.chained_assignment", None):
        sys.exit(main(*sys.argv[1:]))
