#! /usr/bin/env python

import os, subprocess
import pysam

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
                  'W': 'W', 'S': 'S', 'M': 'K', 'K': 'M', 'R': 'Y', 'Y': 'R',
                  'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D', 'N': 'N'}

    return("".join(complement.get(base, base) for base in reversed(seq)))

def extract_read_id(alu_file, bam_file, output_file):

    bam_h = pysam.AlignmentFile(bam_file)
    rname2key = {}
    with open(alu_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            tchr, tstart, tend, tlabel, _, tstrand, _, _, _, ttsd = line.rstrip('\n').split('\t') 
            key = ','.join([tchr, tend, tlabel, tstrand, ttsd])
            print(key)
            for read in bam_h.fetch(tchr, max(int(tend) - 500, 0), int(tend) + 500):

                if read.is_unmapped or read.is_duplicate or read.is_supplementary or read.is_secondary: continue
                if read.mapq < 40: continue
                left_clipping, right_clipping = 0, 0

                # get the clipping size in the both side
                if len(read.cigar) > 1:
                    left_clipping = (read.cigar[0][1] if read.cigar[0][0] in [4, 5] else 0)
                    right_clipping = (read.cigar[len(read.cigar) - 1][1] if read.cigar[len(read.cigar) - 1][0] in [4, 5] else 0)

                if left_clipping >= 8 or right_clipping >= 8 or not read.is_proper_pair: 
                    if read.qname not in rname2key: rname2key[read.qname] = []
                    rname2key[read.qname].append(key)
 

    # remove duplicated keys
    for rname in rname2key:
        keys = list(set(rname2key[rname]))
        rname2key[rname] = keys

    with open(output_file + ".unsorted", 'w') as hout:
        for read in bam_h.fetch():
            if read.is_duplicate or read.is_supplementary or read.is_secondary: continue
            
            if read.qname in rname2key:
                read_seq, read_qual = read.query_sequence, read.qual
                if read.is_reverse: 
                    read_seq = read_seq = reverse_complement(read_seq)
                    read_qual = read_qual[::-1]
                read_qname = read.qname + '/' + ('1' if read.is_read1 else '2')
                for key in rname2key[read.qname]:
                    print(key + '\t' + read_qname + '\t' + read_seq + '\t' + read_qual, file = hout)

    bam_h.close()

    with open(output_file, 'w') as hout:
        subprocess.check_call(["sort", "-k1,1", output_file + ".unsorted"], stdout = hout)
    os.remove(output_file + ".unsorted")



if __name__ == "__main__":

    import sys
    alu_file = sys.argv[1]
    bam_file = sys.argv[2]
    output_file = sys.argv[3]

    extract_read_id(alu_file, bam_file, output_file)

