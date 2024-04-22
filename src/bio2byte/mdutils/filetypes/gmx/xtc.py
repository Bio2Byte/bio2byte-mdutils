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

    XTC_FRAME_FORMAT = (
        ">"
        "l"             # a) Magic number (4) == 1995
        "l"             # b) Atoms (4)
        "l"             # c) Step (4)
        "f"             # d) Time (4)
        "fffffffff"     # e-m) Box (9*4)
        "l"             # n) Atoms (4), checked to be equal to b)
        "f"             # o) Precision (4)
        "llllll"        # p-u) Extent (6*4), MIN: x,y,z, MAX: x,y,z
        "l"             # v) Size in bytes (4)
        # ... Compressed coordinates
    )

    def __init__(self, filename: str):
        self.filename = filename
        self._fhandle = open(filename, "rb")
        self._check()
        self._pointers = list(self._iter_frame_pointers())
        self.goto(0)
        (
            _,              # Magic number
            self.natoms,    # Number of atoms in the system
            self.step0,     # First step in the XTC file
            self.time0      # First time in the XTC file [ps]
        ) = struct.unpack(self.XTC_FRAME_FORMAT[:5], self.read(16))
        self.goto(-1)
        (
            _,              # Magic number
            self.natoms,    # Number of atoms in the system
            self.step1,     # Last step in the XTC file
            self.time1      # Last time in the XTC file [ps]
        ) = struct.unpack(self.XTC_FRAME_FORMAT[:5], self.read(16))
        self.goto(0)        # Reset pointer

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
    
    def seek(self, offset: int, whence: int=0):
        """Goto position in XTC file"""
        self._fhandle.seek(offset, whence)
        return self._fhandle.tell()

    def read(self, n_bytes=None):
        """Read n_bytes starting from the current file location"""
        return self._fhandle.read(n_bytes)
    
    def tell(self):
        """Returns the current position in the XTC file"""
        return self._fhandle.tell()

    def _iter_frame_pointers(self, bytebuffer=10007):
        """Iterates over the file and locates the positions of all frame headers"""
        pos = 0
        self.seek(pos)
        tag = self.read(8)  # Fixed tag: Magic number + Number of atoms
        
        self.seek(pos)
        buffer = self.read(bytebuffer)        
        while len(buffer) >= 56:
            # Search tag in buffer
            idx = buffer.find(tag)
            # If tag was found at end of buffer, extend buffer to include n)
            if len(buffer) < idx+56:
                buffer += self.read(56)
            # If tag was found and b) == n), yield the position and move on
            if idx > -1 and buffer[idx+52:idx+56] == tag[4:]:
                pos += idx
                yield pos
                buffer, pos = buffer[idx+56:], pos + 56
            # If tag was not found, move on
            elif len(buffer) >= 56:
                buffer, pos = buffer[-55:], pos + len(buffer) - 55
            # Extend reading buffer if needed
            if len(buffer) < 56:
                buffer += self.read(bytebuffer)

    def _check(self):
        """Checks if file adheres to XTC format"""
        last_position = self.tell()
        buffer = self.read(struct.calcsize(self.XTC_FRAME_FORMAT))
        header = struct.unpack(self.XTC_FRAME_FORMAT, buffer)
        # Check magic number
        if header[0] != 1995:
            raise XTCParsingError("Unexpected magic number encountered: {0:d}".format(
                header[0]))
        # Check number of atoms, match between b) and n)
        elif header[1] != header[13]:
            raise XTCParsingError("Atom numbers do not match: {0:d} != {1:d}".format(
                header[1], header[13]))
        self.seek(last_position)


def main(args):
    sys.stderr.write("#File\tNumOfAtoms\tNumOfFrames\n")
    for fname in args:
        xtc = GromacsXTC(fname)
        sys.stdout.write("{0}\t{1:d}\t{2:d}\n".format(fname, xtc.natoms, xtc.nframes))

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))