#! /usr/bin/env python

import sys

input_file = sys.argv[1]

with open(input_file, 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        F[1] = str(int(F[1]) - 30)
        F[2] = str(int(F[2]) + 30)
        print('\t'.join(F))

