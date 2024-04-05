#!/usr/bin/env python
import sys
import argparse
import operator as op
import textwrap

import MDAnalysis as mda
import numpy as np
import pandas as pd
import tqdm

def parse_arguments(arguments):
    """ Parse command line arguments passed by th user """

    aparser = argparse.ArgumentParser(description=textwrap.dedent("""\
        Script to analyze ionic contacts in a md simulation.

        Note: Ionic residues are identified automatically by their net 
        charge being != 0."""))

    # Input/output
    aparser.add_argument("--structure", "-s", 
        help=("Structure containing topology information for the trajectory."
              "[tpr, pdb, gro, top]"))
    aparser.add_argument("--trajectory", "-f", nargs="+",
        help=("MD trajectory containing the coordinates to analyze. "
              "[xtc, trr, nc, mdcrd]"))
    aparser.add_argument("--outfile", "-o", 
        help=("(Optional) File to write the output to. "
              "Default: ionic_interactions.csv"))
    
    # Options
    aparser.add_argument("--atom-selection1", "-a",
        help=("(Optional) Atomic selection(s) between which the ionic contacts"
              "are calculated. Default: 'Protein'"))
    aparser.add_argument("--atom-selection2", "-b",
        help=("(Optional) Atomic selection(s) between which the ionic contacts"
              "are calculated. Default: 'Protein'"))
    aparser.add_argument("--distance-cutoff", "-c", type=float, default=5.0,
        help=("(Optional) Maximum distance in Angstrom between two atoms, to still "
              "count as contacts. Default: 5.0"))
    
    args = aparser.parse_args(arguments)
    return args


def find_ionicCenterOfResidue(residue, is_positive=True):
    """ Find "ionic atom(s)" in the residue by finding the heavy atom,
    which has the highest/lowest charge in the atom """
    
    if is_positive:
        comp_q = op.gt
    else:
        comp_q = op.lt
    
    max_q = 0
    atom_list  = []
    
    for hatom in residue.atoms:
        
        if hatom.element == "H" or hatom.type == "C":
            continue
            
        # Get charge of heavy atom + associated hydrogens
        q = hatom.charge + sum((a.charge for a in hatom.bonded_atoms if a.element == "H"))
        if comp_q(q, max_q):
            atom_list = [hatom]
            max_q = q
        elif q == max_q:
            atom_list.append(hatom)
    
    return mda.AtomGroup(atom_list)


def find_ionicAtoms(mddata, selection="*"):
    """ Find ionic atoms in the MD structure """
    
    cations = mda.AtomGroup([], mddata)
    anions = mda.AtomGroup([], mddata)
    
    sele = mddata.select_atoms(selection)
    for res in sele.residues:
        if res.charge > 1e-6: # To account for floating point inaccuracies
            cations = cations.concatenate(find_ionicCenterOfResidue(res, True))
        elif res.charge < -1e-6:
            anions = anions.concatenate(find_ionicCenterOfResidue(res, False))
        else:
            pass        
    return cations, anions


def contacts_within_cutoff(mddata, sele0, sele1, distance_cutoff=5.):
    
    idx0 = sele0.indices
    idx1 = sele1.indices
    contacts = []
    
    for coords in tqdm.tqdm(mddata.trajectory):
        t = coords.time
        crd0 = np.repeat(coords[idx0], idx1.size, axis=0) # slow moving idx
        crd1 = np.tile(coords[idx1].T, idx0.size).T       # fast moving idx
        
        dist = np.sqrt(np.sum((crd1 - crd0) ** 2, axis=1))
        for i in np.where(dist < distance_cutoff)[0]:
            d = dist[i]
            i0 = idx0[i // idx1.size]
            i1 = idx1[i %  idx1.size]
            atm0 = mddata.atoms[i0]
            atm1 = mddata.atoms[i1]
            contacts.append([
                t, atm0.resid, 1+atm0.index, atm1.resid, 1+atm1.index, d
            ])
            
    df = pd.DataFrame(contacts, 
        columns=["time", "cation_resid", "cation_atomid", "anion_resid", "anion_atomid", "distance[A]"],
    )
    return df


def main():
    # Parse commandline arguments
    args = parse_arguments(sys.argv[1:])

    # Load Universe
    MD = mda.Universe(args.structure, args.trajectory)

    # Analyze ionic contacts
    catA, aniA = find_ionicAtoms(MD, args.atom_selection1)
    catB, aniB = find_ionicAtoms(MD, args.atom_selection2)

    contacts = pd.DataFrame()
    for cion, aion in [(catA, aniB), (catB, aniA)]:
        contacts = contacts.append(
            contacts_within_cutoff(MD, cion, aion, 
                distance_cutoff=args.distance_cutoff)
        )

    contacts.to_csv(args.outfile)
    return 0


if __name__ == "__main__":
    sys.exit(main())