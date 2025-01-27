# Matthew Holman            12-15-2019
# Python pound define changer
#
# USE: python changeVals.py ../<directory>/const.h <MATRIX_DIM> <NUM_THREADS>

import re
import sys

f=open(sys.argv[1], "r+")

if f.mode == 'r+':
    contents = f.read()
    contents = re.sub(r"(#define)\s*(MATRIX_DIM)\s*(\d*)", r"#define MATRIX_DIM " + str(sys.argv[2]), contents)
    f.seek(0)
    f.write(contents)
    f.truncate()

f.close()
