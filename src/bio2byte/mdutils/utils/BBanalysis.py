import weakref

import numpy as np
import pandas as pd
from scipy.interpolate import interpn

 

class BBHistogram:
    """ Makes histograms from data in parent BBDihedrals object """
    
    def __init__(self, parent, zero_to_nan=False):
        self.parent         = weakref.ref(parent)
        self.zero_to_nan    = zero_to_nan

    @property
    def histogram_bins(self):
        return self.parent().histogram_bins
        
    def __getitem__(self, resid):
        try:
            phi,psi     =   self.parent()._get_resid(resid)
        except ValueError:
            raise IndexError("residue %d not found." % resid)
        hist        =   self.make_histogram2d(phi,psi)
        return hist
        
    def make_histogram2d(self, xdata, ydata):
        hist,bx,by  =   np.histogram2d(xdata, ydata,
                            bins=self.histogram_bins, density=False)
        if self.zero_to_nan:
            hist[hist == 0] = np.nan
        hist /= xdata.size
        return hist


class BBColoredScatter:
    """ Makes density-colored scatter plots from data in parent
    BBDihedrals object """
    
    def __init__(self, parent, sort=True, wrapped=True, interpolation_method='linear'):
        self.parent = weakref.ref(parent)
        self.sorted = sort
        self.wrapped = wrapped
        self.interpolation_method = interpolation_method

    @property
    def histogram_bins(self):
        return self.parent().histogram_bins
        
    def __getitem__(self, resid):
        phi,psi     =   self.parent()._get_resid(resid)
        dens,xe,ye  =   self._make_histogram2d(phi, psi)
        xyz         =   self._assign_density(phi, psi, dens,
                            xedges=xe, yedges=ye)
        return xyz
        
    def _make_histogram2d(self, xdata, ydata):
        hist,bx,by  =   np.histogram2d(xdata, ydata,
                            bins=self.histogram_bins, density=False)
        hist /= xdata.size
        return hist,bx,by
    
    def _assign_density(self, xdata, ydata, densities, xedges, yedges):
        """
        Assigns an interpolated density to every (x,y) pair
        """
        x_arr = np.array(xdata)
        y_arr = np.array(ydata)
        dens_padded, xe, ye = self._pad_histogram2d(
            densities, xedges, yedges)
        z_arr  = interpn(
            (.5*(xe[:-1]+xe[1:]), .5*(ye[:-1]+ye[1:]),), 
            dens_padded, np.column_stack([x_arr, y_arr]), 
            method=self.interpolation_method, bounds_error=True
        )
        if (z_arr < 0).any(): 
            raise ValueError("Negative density.")
        if np.isnan(z_arr).any():
            raise ValueError("Unexpected NaNs found.")
        if self.sorted:
            idx = z_arr.argsort()
            x_arr, y_arr, z_arr = x_arr[idx], y_arr[idx], z_arr[idx]
        return np.vstack([x_arr, y_arr, z_arr])
    
    def _pad_histogram2d(self, hist, xe, ye, pad_width=1):
        """ Pads 2d histogram (and bin edges) for use with 
        scipy.interpolate.interpn """
        if self.wrapped:
            _hist = np.pad(hist, pad_width=pad_width, mode='wrap')
        else:
            _hist = np.pad(hist, pad_width=pad_width, mode='edge')
        _x    = np.pad(xe, pad_width=pad_width, 
                       mode='reflect', reflect_type="odd")
        _y    = np.pad(ye, pad_width=pad_width, 
                       mode='reflect', reflect_type="odd")    
        return _hist, _x, _y
        

class BBDihedrals:
    """ Stores information of backbone dihedrals """
    
    def __init__(self,  filename, histogram_bins=np.arange(-180,181,9)):
        self._srcfile   =   filename
        self._data      =   pd.read_csv(filename,
                                        index_col=0, header=0)
        
        self.n_values   =   self._data.shape[0]
        self.n_residues =   sum(
                1 for col in self._data.columns 
                if col.endswith("phi") or col.endswith("psi")
            ) // 2
        self.histogram_bins = histogram_bins

        self.scatter = BBColoredScatter(self)
        self.hist    = BBHistogram(self)

        
    
    def __getitem__(self, resid):
        return self._get_resid(resid)
        
    def _get_resid(self, resid):
        """ get data for a given residue number """
        
        columns = [
            col for col in self._data.columns 
            if int(col.split(':')[1]) == resid
        ]
        dih1,dih2,*chi = columns
        if dih1.endswith("phi") and dih2.endswith("psi"):
            phidata = np.array(self._data[dih1])
            psidata = np.array(self._data[dih2])
        elif dih1.endswith("psi") and dih2.endswith("phi"):
            phidata = np.array(self._data[dih2])
            psidata = np.array(self._data[dih1])
        else:
            raise ValueError("data for phi/psi not available in" + str(columns))
        
        return phidata,psidata
    
    def iterresids(self):
        columns = [
            col for col in self._data.columns 
            if col.endswith("phi") or col.endswith("psi")
        ]
        columns = sorted(columns, key=lambda k: tuple(k.split(':')[1:]))
        
        phicol=tuple()
        for col in columns:
            csplit = tuple(col.split(':')[:2])
            if col.endswith("phi"):
                phicol = csplit
            elif col.endswith("psi") and csplit == phicol:
                yield int(phicol[1])   


def boltzmann_energy(probArray, T=298, relative=True):
    arr = probArray.copy()
    arr[arr == 0] = np.nan
    #if np.nansum(probArray) != 1:
    #    probArray /= np.nansum(probArray)
    #R  = 1.9872036e-3   # kcal/(mol K)
    R  = 8.31446261815324e-3 # kJ/(mol K)
    N  = 6.02214076e23  # 1/mol
    dG = - (T*R) * np.log(arr)
    if relative:
        return dG - np.nanmin(dG)
    else:
        return dG


def stack_uneven(arrays, fill_value=np.nan):
    '''
    Fits arrays into a single numpy array, even if they are
    different sizes. `fill_value` is the default value.

    Args:
            arrays: list of np arrays of various sizes
                (must be same rank, but not necessarily same size)
            fill_value (float, optional):

    Returns:
            np.ndarray
    '''
    sizes = [a.shape for a in arrays]
    max_sizes = np.max(list(zip(*sizes)), -1)
    # The resultant array has stacked on the first dimension
    result = np.full((len(arrays),) + tuple(max_sizes), fill_value)
    for i, a in enumerate(arrays):
      # The shape of this array `a`, turned into slices
      slices = tuple(slice(0,s) for s in sizes[i])
      # Overwrite a block slice of `result` with this array `a`
      result[i][slices] = a
    return result
