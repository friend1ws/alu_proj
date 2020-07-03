#! /usr/bin/env python

import sys, gzip
import pysam

mei_file = sys.argv[1]

vcf_file = pysam.VariantFile(mei_file)
for record in vcf_file.fetch():

    lchr = 'chr' + record.contig   
    lstart, lend = str(record.pos - 1), str(record.pos)
    lid = record.id
    ltype, lastart, laend, lstrand = record.info["MEINFO"]
    llen = record.info["MEINFO"][2]
    ltsd = record.info["TSD"]

    lgt = record.samples["NA12878"]["GT"]
    if lgt[0] == 0 and lgt[1] == 0: continue

    # label = ','.join([lchr, lstart, lend, lstrand, lid])    
    label = lid
    laf = float(record.info["AF"][0])
    if laf < 0.001: continue

    print('\t'.join([lchr, lstart, lend, label, str(laf), lstrand, ltype, str(lastart), str(laend), ltsd]))


    
