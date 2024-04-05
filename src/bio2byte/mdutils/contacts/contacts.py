from multiprocessing import Pool

import MDAnalysis as mda 
import numpy as np
import tqdm


class ContactAnalyzer():

    def __init__(self, mddata: mda.Universe, selection0: mda.AtomGroup, selection1: mda.AtomGroup, distance_cutoff: float=6., by_residue: bool=False, nprocs: int=1):
        self._u = mddata
        self._selection0 = selection0
        self._selection1 = selection1
        self._cutoff = distance_cutoff
        self._by_residue = by_residue
        self._nprocs = nprocs

    def _make_sparse_pairs(self):
        lst0 = list(self._selection0.indices)
        lst1 = list(self._selection1.indices)
        i, j =  [], []
        for x in lst0:
            try:
                lst1.remove(x)
            except ValueError: 
                pass
            i += [x] * len(lst1)
            j += lst1
        idcs = np.stack([i,j], axis=1)
        self._sparse_pairs = idcs

    def _contacts_per_frame(self, frame_idx):
        contacts = []
        ts = self._u.trajectory[frame_idx]
        t = ts.time
        dist = np.sqrt(np.sum(
            (ts[self._sparse_pairs[:,0]] - ts[self._sparse_pairs[:,1]]) ** 2, 
            axis=1))

        for i in np.where(dist < self._cutoff)[0]:
            d = dist[i]
            i0, i1 = self._sparse_pairs[i]
            atm0, atm1 = self._u.atoms[i0], self._u.atoms[i1]
            contacts.append([
                t, atm0.resid, 1+atm0.index, atm1.resid, 1+atm1.index, d
            ])

        return contacts

    def _contacts_per_frame_byres(self, frame_idx):
        
        ts = self._u.trajectory[frame_idx]
        t = ts.time
        dist = np.sqrt(np.sum(
            (ts[self._sparse_pairs[:,0]] - ts[self._sparse_pairs[:,1]]) ** 2, 
            axis=1))
        
        contacts = []
        for i in np.where(dist < self._cutoff)[0]:
            d = dist[i]
            i0, i1 = self._sparse_pairs[i]
            r0, r1 = self._u.atoms[i0].resid, self._u.atoms[i1].resid

            try:
                assert r1 != _r1 or r0 != _r0
                contacts.append([t, r0, i0 + 1, r1, i1 + 1, d])
            except NameError:
                contacts.append([t, r0, i0 + 1, r1, i1 + 1, d])
            except AssertionError:
                if d < _d:
                    contacts[-1] = [t, r0, i0 + 1, r1, i1 + 1, d]
                else:
                    continue
            finally:
                _r0, _r1, _d = r0, r1, d

        return contacts


    def run(self):
        """ """
        if self._by_residue:
            func = self._contacts_per_frame_byres
        else:
            func = self._contacts_per_frame

        # Create a sparse list of atomic pairs
        self._make_sparse_pairs()

        # Find contacts
        contacts = []
        frames = np.arange(self._u.trajectory.n_frames)
        with Pool(self._nprocs) as P:
            for cnt in tqdm.tqdm(P.imap(func, frames, chunksize=100)):
                contacts += cnt
        
        self.result = contacts
        return self.result