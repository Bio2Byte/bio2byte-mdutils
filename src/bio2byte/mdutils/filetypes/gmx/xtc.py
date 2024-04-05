#!/usr/bin/env python

"""
Module for parsing GROMACS XTC files in Python.

This module is loosly based on the work of Tsjerk: 
    https://github.com/Tsjerk/xtc/blob/master/xtc.py

It uses a fixed byte format to find the individual frames in
an XTC file:
     0. Magic number (4)    0 -  4 l
     1. Atoms (4)           4 -  8 l
     2. Step  (4)           8 - 12 l 
     3. Time (4)           12 - 16 f
     4. Box (9*4)          16 - 52 fffffffff
    13. Atoms (4)          52 - 56 l           (checked to be equal to b.)
    14. Precision (4)      56 - 60 f
    15. Extent (6*4)       60 - 84 llllll      (MIN: x,y,z, MAX: x,y,z)
    21. smallidx (4)       84 - 88 l
    22. Size in bytes (4)  88 - 92 l 
    23. Coordinates (compressed)
"""

import os
import sys
import struct


class GromacsXTC:

    XTC_FRAME_FORMAT = ">3l10flf8l"

    def __init__(self, filename):
        self.filename = filename
        self.__fhandle = open(filename, "rb")
        n, m, s0, s1, t0, t1 = self._get_trjattr()
        self.__nframes = n
        self.__natoms = m
        self.__step0 = s0
        self.__step1 = s1
        self.__time0 = t0
        self.__time1 = t1

    @property
    def nframes(self):
        ''' Returns total number of frames in the xtc file '''
        return self.__nframes

    @property
    def natoms(self):
        ''' Returns number of atoms in the xtc file '''
        return self.__natoms

    @property
    def step0(self):
        ''' Returns the step number of the first frame '''
        return self.__step0

    @property
    def step1(self):
        ''' Returns the step number of the last frame '''
        return self.__step1

    @property
    def time0(self):
        ''' Returns the time index of the first frame '''
        return self.__time0

    @property
    def time1(self):
        ''' Returns the time index of the last frame '''
        return self.__time1

    def __del__(self):
        self.__fhandle.close()

    def _get_trjattr(self):
        m0, s0, t0 = struct.unpack(">llf", self.read(12, 4))
        n = 0
        for i in self._iter_frame_indices(): n += 1
        m1, s1, t1 = struct.unpack(">llf", self.read(12, i+4))
        return n, m1, s0, s1, t0, t1

    def read(self, n_bytes, pointer=None):
        if isinstance(pointer, int):
            self.__fhandle.seek(pointer)
        return self.__fhandle.read(n_bytes)

    def _iter_frame_indices(self, bytebuffer=10007):
        self.__fhandle.seek(0)
        tag = self.__fhandle.read(8)    # Fixed tag (MagicNumber + NumberOfAtoms)
        self.__fhandle.seek(0)
        buffer = self.__fhandle.read(bytebuffer)
        pos = 0
        while len(buffer) >= 56:
            idx = buffer.find(tag)      # Find tag in buffer
            try:
                # Assert that there is a match
                assert idx > -1

                # Verify match by checking if NumberOfAtoms is repeated at the right position
                if len(buffer) < idx+56:
                    buffer += self.__fhandle.read(56)
                assert buffer[idx+52:idx+56] == tag[4:]
                
                # Match found: yield position and move on
                pos += idx
                yield pos
                buffer = buffer[idx+56:]
                pos += 56

            except AssertionError:
                # No match found: move on
                if len(buffer) >= 56:
                    pos += len(buffer) - 55
                    buffer = buffer[-55:]
            if len(buffer) < 56:
                buffer += self.__fhandle.read(bytebuffer)

def main(*xtcfiles):
    """Displays the number of atoms and number of frames of an XTC file."""
    for xtcfile in xtcfiles:
        if os.path.isfile(xtcfile):
            xtc = GromacsXTC(xtcfile)
            sys.stdout.write("%s\t%d\t%d\n" % (xtcfile, xtc.natoms, xtc.nframes))
        else:
            sys.stderr.write("%s\tERROR: File not found!\n" % xtcfile)

if __name__ == "__main__":
    main(*sys.argv[1:])
    
        
        
    

