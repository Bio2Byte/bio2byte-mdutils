from functools import cached_property
import os
import subprocess

import MDAnalysis as mda
from pytrr import GroTrrReader

from ..utils.shabc import ShellInterface


env_GMX_MAXBACKUP = os.getenv("GMX_MAXBACKUP")
os.environ["GMX_MAXBACKUP"] = "-1"


class GmxDihedralPCA(ShellInterface):
    
    def __init__(self, args, **kwargs):
        self._struct   = args.struct
        self._traj     = args.traj

        if args.readn and not os.path.isfile(args.index or ""):
            raise ValueError(f"Cannot read dihedral index file; file not found: {args.index}")
        if args.readv and not os.path.isfile(args.eigenvec or ""):
            raise ValueError(f"Cannot read eigenvectors; file not found: {args.eigenvec}")
        self._readn = args.readn
        self._readv = args.readv
        self.first_proj = args.first
        self.last_proj = args.last
        
        super().__init__(**kwargs)
        self._dihe_ndx     = self.tmpfile(args.index, suffix=".ndx")
        self._angdist_xvg  = self.tmpfile(args.angdist,  suffix=".xvg")
        self._dihe_trr     = self.tmpfile(args.dihedral_trr, suffix=".trr")
        self._dihe_gro     = self.tmpfile(args.dihedral_gro, suffix=".gro")
        self._eigenval_xvg = self.tmpfile(args.eigenval, suffix=".xvg")
        self._eigenvec_trr = self.tmpfile(args.eigenvec, suffix=".trr")
        self._proj_xvg     = self.tmpfile(args.proj, suffix=".xvg")
        
    @property
    def structure(self):
        return self._struct
    
    @property
    def trajectory(self):
        return self._traj
    
    @cached_property
    def universe(self):
        return mda.Universe(self.structure, self.trajectory)
    
    def run(self):
        # Get the indices for the phi/psi dihedrals
        if not self._readn: self._make_dihedral_ndx()
        # Write phi and psi angles out as sin/cos pairs in a trr-file
        self._make_dihedral_trr()        
        # Make a topology file for the phi/psi dihedrals
        self._make_dihedral_gro()
        # Diagonalize covariance matrix
        if not self._readv: self._covar()
        # Eigenvector analysis
        self._anaeig()
    
    def _make_dihedral_ndx(self):
        phi_indices = []
        psi_indices = []
        
        n_dihedrals = 0
        for res in self.universe.residues:
            phi_atoms = res.phi_selection()
            psi_atoms = res.psi_selection()
            if phi_atoms and psi_atoms:
                phi_indices.append(list(map(str, phi_atoms.indices + 1)))
                psi_indices.append(list(map(str, psi_atoms.indices + 1)))
                n_dihedrals += 2
        self.n_dihedrals = n_dihedrals
        
        with open(self._dihe_ndx, "w") as fhandle:

            fhandle.write("[ Phi_Psi_dihedrals ]\n")
            for phi, psi in zip(phi_indices, psi_indices):
                fhandle.write(" ".join(phi + psi) + "\n")

            fhandle.write("[ Phi_dihedrals ]\n")
            for phi in phi_indices:
                fhandle.write(" ".join(phi) + "\n")

            fhandle.write("[ Psi_dihedrals ]\n")
            for psi in psi_indices:
                fhandle.write(" ".join(psi) + "\n")
        
        return self

    def _make_dihedral_trr(self):
        cmd = ["gmx", "angle", "-f", self.trajectory, "-n", self._dihe_ndx, 
               "-od", self._angdist_xvg, "-or", self._dihe_trr,
               "-type", "dihedral"]
        subprocess.run(cmd, input=b"0\n")
        return self
    
    def _make_dihedral_gro(self):
        # For each dihedral cos and sin are stored, while an atom can store up to three values
        with GroTrrReader(self._dihe_trr) as trr:
            header, _ = trr.read_frame(read_data=False)
        n_atoms = header["natoms"]
             
        with open(self._dihe_gro, "w") as fhandle:
            fhandle.write(f"Structure file for dihedral-PCA containing {self.n_dihedrals} dihedrals\n")
            fhandle.write(f"{n_atoms:6d}\n")
            for idx in range(n_atoms):
                resid = (idx // 4) + 1
                atmid = min(99999, idx + 1)
                resname = "ALA"
                atmname = "CA"
                fhandle.write(f"{resid:5d}{resname:>5}{atmname:>5}{atmid:5d}{'   0.000'*3}{'  0.0000'*3}\n")
        
        return self
    
    def _covar(self):
        cmd = ["gmx", "covar", 
               "-s", self._dihe_gro, "-f", self._dihe_trr, 
               "-o", self._eigenval_xvg, "-v", self._eigenvec_trr,
               "-av", self.tmpfile(None, suffix=".pdb"),     # Absolutely useless average structure
               "-l", self.tmpfile(None, suffix=".log"),      # Absolutely useless logfile
               "-xvg", "none", "-nofit"]
        subprocess.run(cmd, input=b"0\n")
        return  self
    
    def _anaeig(self):
        cmd = ["gmx", "anaeig", "-v", self._eigenvec_trr, "-f", self._dihe_trr, "-s", self._dihe_gro, 
               "-proj", self._proj_xvg, "-xvg", "none", "-first", self.first_proj, "-last", self.last_proj]
        subprocess.run(cmd, input=b"0\n0\n")
        return self
