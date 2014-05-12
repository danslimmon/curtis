#!/usr/bin/python

import sys
import math
import struct
import binascii

import bitarray
from scipy import zeros, sin, pi
from scipy.fftpack import dct, idct


class CompressedDataset:
    def __init__(self):
        self.bits = bitarray.bitarray()

    def append_block(self, new_bits):
        """Extends the given bitarray onto our data."""
        self.bits.extend(new_bits)

    def tofile(self, f):
        """Writes the compressed data to the given file."""
        self.bits.tofile(f)

    def fromfile(self, f):
        """Reads compressed data from the given file."""
        self.bits.fromfile(f)

    def decompress(self, qmat, block_size,
                   point_count, interval, min_ts):
        """Decompresses the bitarray and returns a Dataset object."""
        coeffs = zeros(block_size)
        values = zeros(point_count)
        ind = 0
        i = 0
        block_num = 0
        while ind < self.bits.length():
            if qmat[i] == 52:
                coeffs[i] = 0.
            elif ind > len(self.bits) - (64 - qmat[i]):
                # File is over
                break
            else:
                v = self.bits[ind:ind + 64 - qmat[i]]
                v.extend(qmat[i] * (False,))
                coeffs[i] = struct.unpack(">d", v.tobytes())[0]
                ind += 64 - qmat[i]

            i += 1

            if i >= block_size:
                values[block_num * block_size:block_num * block_size + block_size] = idct(coeffs, norm='ortho')
                i = 0
                block_num += 1
                # We pad out to a full byte at the end of the block
                if ind % 8 != 0:
                    ind += 8 - (ind % 8)
                if i > 1:
                    raise SystemExit

        return Dataset(point_count, interval, min_ts, values)


class Dataset:
    def __init__(self, point_count, interval, min_ts, values):
        self.point_count = point_count
        self.interval = interval
        self.min_ts = min_ts
        self.values = values

    def compress(self, qmat, block_size):
        """Runs compression on the data and returns an array, returning a CompressedDataset instance."""
        compressed = CompressedDataset()
        for i in range(0, self.point_count, block_size):
            if i + block_size > self.point_count:
                break
            d = dct(self.values[i:i+block_size], norm='ortho')
            compressed.append_block(self._quantize(d, qmat))
        return compressed

    def _quantize(self, block, qmat):
        """Quantizes each value in the frequency-domain block according to the quantization matrix.
        
           Returns a bigint representation of the quantized block."""
        bits = bitarray.bitarray()
        for i in range(len(block)):
            if qmat[i] == 52:
                continue
            # Extend `bits` with the latest float, encoded as a double
            bits.frombytes(struct.pack('>d', block[i]))
            for i in range(qmat[i]):
                bits.pop()
        bits.fill()
        return bits

def gen_qmat(n):
    """Generates a constant quantization matrix to be used on the frequency components.

       It's an array of integers, each representing the number of bits to mask off of the
       corresponding floating-point value's mantissa.
 
       Values in the quantization matrix range from 0 to 52, the size of a double-precision
       float's significand. A value of 52 indicates dropping the value entirely."""
    m = zeros(n, dtype=int)
    for i in range(n):
        m[i] = min(52, int(53 * (float(i*2)/float(n))**3))
    return tuple(m)

def read_whisperdump(path):
    """Turns a whisper dump file into a Dataset object."""
    f = open(path)
    min_ts = 1 << 63
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


QMAT = gen_qmat(2048)
print "Reading dump"
ds = read_whisperdump(sys.argv[1])
print "Compressing"
cds = ds.compress(QMAT, 2048)
print "Writing"
f = open('/tmp/tru.dat', 'wb')
cds.tofile(f)
f.close()

f = open('/tmp/tru.dat', 'rb')
print "Decompressing"
cds = CompressedDataset()
cds.fromfile(f)
ds2 = cds.decompress(QMAT, 2048, ds.point_count, ds.interval, ds.min_ts)
f.close()

f = open("/tmp/curtis.csv", "w")
f.write("pre,post\n")
for v1, v2 in zip(ds.values, ds2.values):
    f.write("{0},{1}\n".format(v1,v2))
