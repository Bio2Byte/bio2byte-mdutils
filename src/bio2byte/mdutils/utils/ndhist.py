from typing import Iterable, Union
import numpy as np


class NDHistogram:
    
    def __init__(self, data, max_bin_limit = 1e5):
        self.data = data
        self._bins = None
        self._ndhists = None
        self.max_bin_limit = max_bin_limit
        
    @property
    def n_replicas(self):
        return self.data.shape[0]
    
    @property
    def bins(self):
        if self._bins is None:
            self.make_nbins()
        return self._bins
    
    @property
    def hists(self):
        if self._ndhists is None:
            self.make_histograms()
        return self._ndhists
    
    @bins.setter
    def bins(self, bins_iterable):
        if len(bins_iterable) != self.data.shape[-1]:
            raise ValueError(f"Data has {self.data.shape[-1]} dimensions; bins cover {len(bins_iterable)} dimensions.")
        else:
            self._bins = bins_iterable
    
    @property
    def ranges(self):
        return np.stack([
            np.min(self.data, axis=(0,1)), 
            np.max(self.data, axis=(0,1))], axis=1)
    
    def iter_histograms(self, first=0, inc_first=0, last=1000, inc_last=1000, **kwargs):
        
        if inc_first == 0 and inc_last == 0:
            raise ValueError("Either 'inc_first' or 'inc_last' must be > 0.")
            
        n = self.data.shape[1] + 1
        for fst, lst in zip(range(first, n, inc_first) if inc_first > 0 else (first for _ in range(last, n, inc_last)),
                            range(last, n, inc_last) if inc_last > 0 else (last for _ in range(first, n, inc_first))):
            self.make_histograms(first=fst, last=lst, **kwargs)
            yield lst-fst, self.hists
    
    def make_histograms(self, *, first=None, last=None, density=False, **kwargs):
        
        # Check max_bin_limit
        n_bins = np.prod([b.size - 1 for b in self.bins])
        if n_bins > self.max_bin_limit:
            raise ValueError(f"Maximum allowed number of bins exceeded: {n_bins:,} bins")
        hists = np.array([
            np.histogramdd(self.data[i,first:last,:], bins=self.bins, density=density, **kwargs)[0]
            for i in range(self.data.shape[0])])
        self._ndhists = hists
        
    def make_nbins(self, nbins: Union[int,Iterable]=10):
        """ Divides the range of each dimension in n bins """
        if isinstance(nbins, int):
            self.bins = [np.linspace(i, j, nbins+1) for i, j in self.ranges]
        elif isinstance(nbins, Iterable):
            self.bins = [np.linspace(i, j, k+1) for (i,j),k in zip(self.ranges, nbins)]
        else:
            raise TypeError(f"n_bins must be of type 'int' or 'Iterable': Got '{nbins}'")
        
    def make_dbins(self, dbins: float = None, *, nbins: int = 10):
        """ Divides the range of each dimension segments of size n """
        if dbins is None and isinstance(nbins, int):
            dbins = np.max(self.ranges[:,1] - self.ranges[:,0]) / nbins
        self.bins = [
            np.linspace(i, j, max(2, int((j - i) / dbins) + 1)) for i, j in self.ranges
        ]
    
    @classmethod
    def from_projectionfile(cls, filename, n_replicas: int, **kwargs):
        ts, prj = np.genfromtxt(filename).T
        ts0, = np.where(ts == ts[0])
        prj = prj.reshape((ts0.size, n_replicas, -1))
        prj = np.moveaxis(prj, 0, -1)
        return cls(prj, **kwargs)
