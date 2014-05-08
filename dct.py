#!/usr/bin/python

import sys
import math
import struct
import binascii

from scipy import zeros, sin, pi
from scipy.fftpack import dct, idct


class Dataset:
    def __init__(self, point_count, interval, min_ts, values):
        self.point_count = point_count
        self.interval = interval
        self.min_ts = min_ts
        self.values = values

    def compress(self, qmat, block_size):
        """Runs compression on the data and returns an array, returning a CompressedDataset instance."""
        compressed = []
        for i in range(0, self.point_count, block_size):
            if i + block_size > len(self.values):
                break
            d = dct(self.values[i:i+block_size], norm='ortho')
            compressed.append(self._quantize(d, qmat))
        return compressed

    def _quantize(self, block, qmat):
        """Quantizes each value in the frequency-domain block according to the quantization matrix.
        
           Returns a bigint representation of the quantized block."""
        bits = 0
        for i in range(len(block)):
            mantissa, exponent = math.frexp(block[i])
            mantissa = floatToBits(mantissa)
            bits <<= 11
            bits |= mantissa
            bits <<= (11 - qmat[i])
            bits |= exponent
        return bits

def floatToBits(x):
    """Returns the given float, reinterpreted as an integer for bitwise math."""
    s = struct.pack('>f', x)
    return struct.unpack('>l', s)[0]

def bitsToFloat(i):
    """Turns the integer (from floatToBits) back into a float."""
    s = struct.pack('>l', i)
    return struct.unpack('>f', s)[0]

def gen_qmat(n):
    """Generates a constant quantization matrix to be used on the frequency components.
    
       It's an array of integers, each representing the number of bits to mask off of the
       corresponding floating-point value's mantissa.
       
       Values in the quantization matrix range from 0 to 52, the size of a double-precision
       float's significand."""
    m = zeros(n, dtype=int)
    for i in range(n):
        m[i] = int(52 * (float(i)/float(n))**3)
    return m

def read_whisperdump(path):
    """Turns a whisper dump file into a Dataset object."""
    f = open(path)
    min_ts = 9999999999999999999999
    in_archive = False
    for line in f:
        if in_archive:
            if line.strip() == "":
                break
            a[i] = float(line.split()[2])
            ts = int(line.split()[1].rstrip(','))
            if ts > 0:
                min_ts = min(ts, min_ts)
            i += 1
            continue
        if line.startswith('  points:'):
            point_count = int(line.split()[1])
            a = zeros(point_count)
        if line.startswith('  seconds per point:'):
            interval = int(line.split()[3])
        if line.startswith('Archive 0 data'):
            i = 0
            in_archive = True
    return Dataset(point_count, interval, min_ts, a)


QMAT = gen_qmat(256)
print "Reading dump"
ds = read_whisperdump(sys.argv[1])
print "Compressing"
cds = ds.compress(QMAT, 256)
print type(cds)
print cds
