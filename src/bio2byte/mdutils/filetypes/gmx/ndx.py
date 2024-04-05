"""
Module for parsing GROMACS index files in python.
"""

import numpy as np

class GromacsNDX:
    """ Class to parse/represent GROMACS index files. """
    
    def __init__(self, filepath, index_correction=-1):
        self._filepath = filepath
        self._indices  = None
        self._group_labels = None
        self.parse_ndx(filepath, index_correction)
        
    def __repr__(self):
        return "<GromacsNDX with %d groups>" % len(self._group_labels)
    
    def __getitem__(self, key):
        if isinstance(key, int):
            return self._indices[self._group_labels[key]]
        elif isinstance(key, str):
            return self._indices[key]
        else:
            raise ValueError("Unrecognized key type: only 'int' or 'str'")
        
    def parse_ndx(self, filepath, index_correction=-1):
        group_labels = []
        indices = {}
        
        with open(filepath, 'r') as fhandle:
            for ln in fhandle:
                if ln.startswith("["):
                    try:
                        group_labels.append(grp)
                        indices[grp] = np.array(idc, dtype=int) + index_correction
                    except NameError:
                        pass
                    grp = ln.split(" ")[1]
                    idc = []
                else:
                    idc.extend(map(int, ln.split()))
            try:
                group_labels.append(grp)
                indices[grp] = np.array(idc, dtype=int) + index_correction
            except NameError:
                pass
        
        self._indices = indices
        self._group_labels = group_labels
