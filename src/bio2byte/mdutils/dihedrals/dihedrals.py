import MDAnalysis as mda
from MDAnalysis.analysis.dihedrals import Dihedral
import numpy as np
import pandas as pd


class DihedralAnalyzer():
    
    def __init__(self, args):
        self.struct =   args.struct
        self.traj   =   args.traj
        self.mask   =   args.mask
        self.unwrap =   args.unwrap

        # Load universe
        self.Unv  =  mda.Universe(self.struct, self.traj,
                                  tpr_resid_from_one=True)

        self.sele =  self.Unv.select_atoms(self.mask)

        if self.unwrap:
            self.Unv.add_transformations(
                mda.transformations.unwrap(self.sele))


    def calculate_dihedrals(self, dihedrals=["phi", "psi", "chi1"]):
        """ Calculates dihedrals """
        columns, atomsele = self._get_dihedrals(dihedrals)

        Dihedrals = Dihedral(atomsele)
        Dihedrals.run()

        df = pd.DataFrame(data = Dihedrals.results.angles,
                          columns = columns)
        return df


    def _get_dihedrals(self, dihedrals=[]):
        def get_phi(residue):
            sele = residue.phi_selection()
            return sele

        def get_psi(residue):
            sele = residue.psi_selection()
            return sele

        def get_chi1(residue):
            for cg in ("SG", "OG", "OG1", "CG", "CG1"):
                sele = residue.chi1_selection(cg_name=cg)
                if sele: break
            return sele

        #def get_chi2(residue):
        #    pass

        get_dihedral = dict(phi  = get_phi,  
                            psi  = get_psi,
                            chi1 = get_chi1)
        columns  = []
        atomsele = []

        for res in self.sele.residues:
            prefix = '{0}:{1:d}'.format(res.resname, res.resid)

            for da in dihedrals:
                sele = get_dihedral[da](res)
                if sele:
                    columns.append('{0}:{1}'.format(prefix, da))
                    atomsele.append(sele)

        return columns, atomsele