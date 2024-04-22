#!/usr/bin/env python

"""
Author:   David Bickel
Date:     22/04/2024

Module for parsing GROMACS XTC files in Python.

This module is loosly based on the work of Tsjerk: 
    https://github.com/Tsjerk/xtc/blob/master/xtc.py
"""

import sys
import struct


class XTCParsingError(ValueError):
    pass


class GromacsXTC:

    # Many thanks to Tsjerk:
    # https://github.com/Tsjerk/xtc/blob/master/xtc.py
    XTC_FRAME_FORMAT = (
        ">"
        "l"             # Magic number (4) == 1995
        "l"             # Atoms (4)
        "l"             # Step (4)
        "f"             # Time (4)
        "fffffffff"     # Box (9*4)
        "l"             # Atoms (4), checked to be equal to above
        "f"             # PRecision (4)
        "lllllll"       # Extent (6*4), MIN: x,y,z, MAX: x,y,z
        "l"             # Size in bytes (4)
        # ... Compressed coordinates
    )

    def __init__(self, filename: str):
        self.filename = filename
        self._fhandle = open(filename, "rb")
        self._pointers = list(self._iter_frame_pointers())
        self.goto(0)
        (
            magicnumber,    # Magic number
            self.natoms,    # Number of atoms in the system
            self.step0,     # First step in the XTC file
            self.time0      # First time in the XTC file [ps]
        ) = struct.unpack(self.XTC_FRAME_FORMAT[:5], self.read(16))
        self.goto(-1)
        (
            magicnumber,    # Magic number
            self.natoms,    # Number of atoms in the system
            self.step1,     # Last step in the XTC file
            self.time1      # Last time in the XTC file [ps]
        ) = struct.unpack(self.XTC_FRAME_FORMAT[:5], self.read(16))
        self.goto(0)        # Reset pointer
        if magicnumber != 1995:
            raise XTCParsingError("Unexpected magic number in XTC file: {0}".format(magicnumber))

    @property
    def nframes(self):
        """Returns total number of frames in the xtc file"""
        return len(self._pointers)

    def __del__(self):
        self._fhandle.close()

    def goto(self, index: int):
        """Goto frame with index in XTC file"""
        self._fhandle.seek(self._pointers[index])
        return self._fhandle.tell()
    
    def seek(self, pos: int):
        """Goto position in XTC file"""
        self._fhandle.seek(pos)
        return self._fhandle.tell()

    def read(self, n_bytes=None):
        """Read n_bytes starting from the current file location"""
        return self._fhandle.read(n_bytes)

    def _iter_frame_pointers(self, bytebuffer=10007):
        pos = 0
        self.seek(pos)
        tag = self.read(8)    # Fixed tag (MagicNumber + NumberOfAtoms)

        self.seek(pos)
        buffer = self.read(bytebuffer)        
        while len(buffer) >= 56:
            idx = buffer.find(tag)      # Find tag in buffer
            try:
                # Assert that there is a match
                assert idx > -1

                # Verify match by checking if NumberOfAtoms is repeated at the right position
                if len(buffer) < idx+56:
                    buffer += self.read(56)
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
                buffer += self.read(bytebuffer)

def main(args):
    sys.stderr.write("#File\tn_atoms\tn_frames\n")
    for fname in args:
        xtc = GromacsXTC(fname)
        sys.stdout.write("{0}\t{1:d}\t{2:d}\n".format(fname, xtc.natoms, xtc.nframes))

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))