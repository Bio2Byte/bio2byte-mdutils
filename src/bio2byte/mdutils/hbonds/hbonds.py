"""
Module with classes that describe hydrogen bonds and mediated hydrogen bonds

* abstract base class  HBondBase
* class HBond
* class MediatedHBond
"""
from abc import ABCMeta, abstractmethod
import functools
from typing import Tuple
import MDAnalysis as mda


class HBondBase(metaclass=ABCMeta):
    """ Base class for all objects representing Hydrogen bonds"""
    
    def __init__(self, universe: mda.Universe, indextuple: Tuple[int]):
        self._u = universe
        self._idxtuple = indextuple

    def __repr__(self):
        return "<HBond {1} frame={0:d}>".format(self.frame, self.to_string())
    
    def __str__(self):
        return self.to_string()
        
    def __eq__(self, other):
        """ Define equality for filtering for duplicates """
        if isinstance(other, self.__class__):
            return self._idxtuple == other._idxtuple
        else:
            return False

    def __lt__(self, other):
        """ Define less-than for sorting """
        if isinstance(other, HBondBase):   
            return self.indextuple < other.indextuple
        else:
            raise TypeError(f"Cannot compare {self.__class__.__name__} with {other.__class__.__name__}")
    
    @property
    def indextuple(self):
        return self._idxtuple
    
    @property
    def frame(self):
        return self.indextuple[0]
    
    @property
    def universe(self):
        return self._u
    
    @property
    @abstractmethod
    def interaction_type(self):
        pass
    
    @abstractmethod
    def to_csv(self):
        pass
    
    @abstractmethod
    def to_string(self):
        pass
    
    def is_backbone(self, atom_index):
        atomname = self.universe.atoms[atom_index].name
        return atomname in ("N", "O", "H", "C")


class HBond(HBondBase):
    """ Class representing a single hydrogen bond """

    def __init__(self, universe: mda.Universe, frame: int, donor_id: int, hydrogen_id: int, acceptor_id: int, distance: float, angle: float, 
                 donor_is_solvent: bool=False, acceptor_is_solvent: bool=False):
        super().__init__(universe, (int(frame), int(donor_id), int(acceptor_id), int(hydrogen_id)))
        self._distance = distance
        self._angle = angle
        self._donor_is_solvent = donor_is_solvent
        self._acceptor_is_solvent = acceptor_is_solvent

    @property
    def distance(self):
        return self._distance

    @property
    def angle(self):
        return self._angle

    @functools.cached_property
    def acceptor(self):
        """ Returns the acceptor atom """
        accid = self._idxtuple[2]
        return self.universe.atoms[accid]

    @functools.cached_property
    def donor(self):
        """ Returns the donor atom """
        donid = self._idxtuple[1]
        return self.universe.atoms[donid]

    @functools.cached_property
    def solvent(self):
        """ Returns the solvent residue """
        if self._donor_is_solvent:
            return self.donor.residue
        elif self._acceptor_is_solvent:
            return self.acceptor.residue
        else:
            raise ValueError(f"No solvent defined for this HBond: {self}")
    
    @functools.cached_property
    def solute(self):
        """ Returns the solute residue """
        if self._donor_is_solvent:
            return self.acceptor.residue
        elif self._acceptor_is_solvent:
            return self.donor.residue
        else:
            raise ValueError(f"No solvent defined for this HBond: {self}")
    
    @functools.cached_property
    def interaction_type(self):
        return "{0}->{1}".format("bb" if self.is_backbone(self.indextuple[1]) else "sc",
                                 "bb" if self.is_backbone(self.indextuple[2]) else "sc")
            
    def to_string(self):
        return "{0}:{1:d}->{2}:{3:d}".format(
            self.donor.resname, self.donor.resid, 
            self.acceptor.resname, self.acceptor.resid)
    
    def to_csv(self):
        return ",".join([
            str(self.frame),
            self.donor.resname, str(self.donor.resid),
            self.acceptor.resname, str(self.acceptor.resid),
            self.interaction_type
        ])


class MediatedHBond(HBondBase):
    """ Class representing hydrogen bonds mediated by a solvent molecule """

    def __init__(self, hbond1: HBond, hbond2: HBond = None):
        self._hbond1, self._hbond2 = sorted([hbond1, hbond2])
        super().__init__(self._hbond1.universe, self._hbond1.indextuple + self._hbond2.indextuple)

    @property
    def distances(self):
        return self._hbond1.distance, self._hbond2.distance

    @property
    def angles(self):
        return self._hbond1.angle, self._hbond2.angle
    
    @property
    def solute1(self):
        return self._hbond1.solute

    @property
    def solvent(self):
        return self._hbond1.solvent
        
    @property
    def solute2(self):
        return self._hbond2.solute
    
    @functools.cached_property
    def interaction_type(self):
        if self._hbond1._donor_is_solvent:
            outstr = "{0}<-aq".format("bb" if self.is_backbone(self.indextuple[2]) else "sc")
        else:
            outstr = "{0}->aq".format("bb" if self.is_backbone(self.indextuple[1]) else "sc")
        if self._hbond2._donor_is_solvent:
            outstr += "->{0}".format("bb" if self.is_backbone(self.indextuple[6]) else "sc")
        else:
            outstr += "<-{0}".format("bb" if self.is_backbone(self.indextuple[5]) else "sc")
        return outstr
    
    def to_string(self):        
        return "{0}:{1:d}{6}{2}:{3:d}{7}{4}:{5:d}".format(
            self.solute1.resname, self.solute1.resid,
            self.solvent.resname, self.solvent.resid,
            self.solute2.resname, self.solute2.resid,
            "<-" if self._hbond1._donor_is_solvent else "->",
            "->" if self._hbond2._donor_is_solvent else "<-")
    
    def to_csv(self):
        return ",".join([
            str(self.frame),
            self.solute1.resname, str(self.solute1.resid),
            self.solute2.resname, str(self.solute2.resid),
            self.interaction_type
        ])