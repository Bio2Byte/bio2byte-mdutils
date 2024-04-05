"""Module for parsing GROMACS XPM files in python."""

import numpy as np

class GromacsXPM:
    
    def __init__(self, filename):
        self.filename = filename
        self.n_rows = None
        self.n_cols = None
        self.n_thresholds = None
        self.thresholds = None
        self.symbols = None
        self.colors = None
        self._update_from_file()

    def _update_from_file(self):
        
        thresholds = []
        symbols = []
        colors = []
        
        with open(self.filename, 'r') as fhandle:
            for ln in fhandle:
                if ln.startswith("static char *gromacs_xpm[] ="):
                    break
            ln = next(fhandle).split('"')[1]
            n_rows, n_cols, n_thresholds, _ = map(int, ln.split())
            for ln in fhandle:
                if not ln.startswith('"'):
                    break
                ln = ln.split('"')
                lns = ln[1].split()
                thresholds.append(float(ln[3]))
                symbols.append(lns[0])
                colors.append(lns[2])
        
        self.n_rows = n_rows
        self.n_cols = n_cols
        self.n_thresholds = n_thresholds
        self.thresholds = np.array(thresholds)
        self.symbols = np.array(symbols, dtype="U1")
        self.colors = np.array(colors, dtype="U7")
        
    def value2symbol(self, value):
        d = np.abs(np.subtract.outer(self.thresholds, value))
        i = np.argmin(d, axis=0)
        return self.symbols[i]