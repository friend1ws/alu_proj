#! /usr/bin/env python

import sys

input_file = sys.argv[1]
output_file1 = sys.argv[2]
output_file2 = sys.argv[3]

hout1 = open(output_file1, 'w')
hout2 = open(output_file2, 'w')

with open(input_file, 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        if F[0] == "chr1,23654652,ALU_umary_ALU_69,+,AACAGATGGGGCATC":

            if F[1].endswith('1'):
                print('@' + F[1] + '\n' + F[2] + '\n' + '+' + '\n' + F[3], file = hout1)
            else:
                print('@' + F[1] + '\n' + F[2] + '\n' + '+' + '\n' + F[3], file = hout2)

hout1.close()
hout2.close()
        

