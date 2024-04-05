import functools

import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis

from .hbonds import HBond, MediatedHBond


class MediatedHBondAnalyzer:
    """ Class to perform hydrogen bond analysis between two selections """
    
    @staticmethod
    def get_donorHAtoms(ag: mda.AtomGroup) -> mda.AtomGroup:
        donH = ag.select_atoms("name H* and bonded name O* N* S*")
        return donH

    @staticmethod
    def get_acceptorAtoms(ag: mda.AtomGroup) -> mda.AtomGroup:
        accO = ag.select_atoms("name O*")
        # Exclude amide Ns (N) and aromatic hydrogens where the lone pair is aromatized (NA)
        accN = ag.select_atoms("name N* and (not type N NA) and (not bonded name O*)" )
        n4list = [N for N in accN.atoms if len(N.bonded_atoms) > 3 or N.charge > -.1]
        if n4list:
            accN -= mda.AtomGroup(n4list)
        return (accO | accN)   
    
    @staticmethod
    def merge_hbonds(hbonds0, hbonds1):
        """ Merges solvent mediated H Bonds """

        med_hbonds = []
        curr_frame, curr_index = -99, 0

        for hb0 in hbonds0:
            
            if hb0.frame > curr_frame:
                curr_frame = hb0.frame
                curr_hbonds1 = []
                for i, hb1 in enumerate(hbonds1[curr_index:]):
                    if hb1.frame == curr_frame:
                        curr_hbonds1.append(hb1)
                    elif hb1.frame > curr_frame:
                        curr_index += i
                        break
                    else:
                        raise ValueError("hbonds1 are not sorted!")                    
            elif hb0.frame < curr_frame:
                raise ValueError("hbonds0 are not sorted!")
            
            for hb1 in curr_hbonds1:
                if hb0.solvent == hb1.solvent and hb0 != hb1:
                    yield MediatedHBond(hb0, hb1)
    
    def __init__(self, universe, sele0, sele1=None, sele_solvent=None, 
                 distance_cutoff=3.0, angle_cutoff=150., update_selection=False):
        self._u = universe
        self._sele0 = sele0
        self._sele1 = sele1 or sele0
        self._solvent = sele_solvent or ""
        self._cutoff_distance = distance_cutoff
        self._cutoff_angle = angle_cutoff
        self._sele_update = update_selection
    
    @property
    def universe(self):
        return self._u

    @property
    def n_frames(self):
        return self._u.trajectory.n_frames
        
    @functools.cached_property
    def donorHs_sele0(self):
        ag = self.__class__.get_donorHAtoms(self._u.select_atoms(self._sele0))
        return ag
    
    @functools.cached_property
    def donorHs_sele1(self):
        ag = self.__class__.get_donorHAtoms(self._u.select_atoms(self._sele1))
        return ag
    
    @functools.cached_property
    def donorHs_solvent(self):
        ag = self.__class__.get_donorHAtoms(self._u.select_atoms(self._solvent or ""))
        return ag
    
    @functools.cached_property
    def acceptors_sele0(self):
        ag = self.__class__.get_acceptorAtoms(self._u.select_atoms(self._sele0))
        return ag
    
    @functools.cached_property
    def acceptors_sele1(self):
        ag = self.__class__.get_acceptorAtoms(self._u.select_atoms(self._sele1))
        return ag
    
    @functools.cached_property
    def acceptors_solvent(self):
        ag = self.__class__.get_acceptorAtoms(self._u.select_atoms(self._solvent or ""))
        return ag
        
    def run(self):
        """ Runs the h-bond analysis """
        
        # Calculate hydrogen-bonds between sele0 and sele1
        direct_hbonds = self._get_direct_hbonds()
        
        # Calculate solvent-mediated hydrogen bonds
        mediated_hbonds = self._get_mediated_hbonds()
        
        all_hbonds = np.concatenate([direct_hbonds, mediated_hbonds], axis=0)
        self.results = all_hbonds[np.argsort(all_hbonds)]
        return self
        
    def _get_direct_hbonds(self):
        """ Runs the h-bond analysis between sele0 and sele1 """
        
        if self._sele0 == self._sele1:
            hbonds = self._get_hbonds(self.donorHs_sele0, self.acceptors_sele1)
        else:
            hbonds = np.concatenate([
                self._get_hbonds(self.donorHs_sele0, self.acceptors_sele1),
                self._get_hbonds(self.donorHs_sele1, self.acceptors_sele0)
                ], axis=0)
        return np.unique(hbonds)

    def _get_mediated_hbonds(self):
        """ Runs the h-bond analysis between sele0 / sele1 and the solvent """
        
        if not self._solvent:
            return np.zeros((0,), dtype="O")
        
        # Calculate hydrogen-bonds between sele0 and the solvent
        slt0_slv = np.unique(np.concatenate([
            self._get_hbonds(self.donorHs_sele0, self.acceptors_solvent, solvent_acceptor=True),
            self._get_hbonds(self.donorHs_solvent, self.acceptors_sele0, solvent_donor=True)], axis=0))
        
        if self._sele0 != self._sele1:
            slt1_slv = np.unique(np.concatenate([
                self._get_hbonds(self.donorHs_sele1, self.acceptors_solvent, solvent_acceptor=True),
                self._get_hbonds(self.donorHs_solvent, self.acceptors_sele1, solvent_donor=True)], axis=0))
        else:
            slt1_slv = slt0_slv
            
        mediated_hbonds = np.array(list(self.__class__.merge_hbonds(slt0_slv, slt1_slv)))
        return np.unique(mediated_hbonds)
    
    def _get_hbonds(self, donorHs, acceptors, solvent_donor=False, solvent_acceptor=False):
        
        if isinstance(donorHs, mda.AtomGroup):
            doh = "index " + " ".join(map(str, donorHs.indices))
        else:
            doh = donorHs
        if isinstance(acceptors, mda.AtomGroup):
            acc = "index " + " ".join(map(str, acceptors.indices))
        else:
            acc = acceptors
        
        Hba = HydrogenBondAnalysis(
            universe      = self.universe,
            donors_sel    = None,
            hydrogens_sel = doh,
            acceptors_sel = acc,
            d_a_cutoff    = self._cutoff_distance,
            d_h_a_angle_cutoff = self._cutoff_angle,
            update_selections = self._sele_update)
        Hba.run()
        
        func = lambda row: HBond(self.universe, *row, solvent_donor, solvent_acceptor)
        return np.apply_along_axis(func, 1, Hba.results.hbonds)