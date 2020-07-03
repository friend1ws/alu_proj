#! /usr/bin/env python

import os, pathlib, subprocess
import pysam

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
                  'W': 'W', 'S': 'S', 'M': 'K', 'K': 'M', 'R': 'Y', 'Y': 'R',
                  'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D', 'N': 'N'}

    return("".join(complement.get(base, base) for base in reversed(seq)))

Alu = "GGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGCGGGCGGATCACGAGGTCAGGAGATCGAGACCATCCTGGCTAACACGGTGAAACCCCGTCTCTACTAAAAATACAAAAAATTAGCCGGGCGTGGTGGCGGGCGCCTGTAGTCCCAGCTACTCGGGAGGCTGAGGCAGGAGAATGGCGTGAACCCGGGAGGCGGAGCTTGCAGTGAGCCGAGATCGCGCCACTGCACTCCAGCCTGGGCGACAGAGCGAGACTCCGTCTCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"

def generate_alu_seq(reference_h, tchr, tpos, tstrand, tsd):

    # variant_seq
    variant_seq = ""
    tseq = reference_h.fetch(region = tchr + ':' + str(int(tpos) - 500 + 1) + '-' + str(int(tpos)))
    if tstrand == '+':
        tseq = tseq + Alu 
    else:
        tseq = tseq + reverse_complement(Alu)
    tseq = tseq + tsd
    tseq = tseq + reference_h.fetch(region = tchr + ':' + str(int(tpos) + 1) + '-' + str(int(tpos) + 500)) 

    return(tseq)


#! /usr/bin/env python

import re, sys, math


def proc_mpileup(input_file, output_file, min_variant_num = 3, min_variant_ratio = 0.25):

    hout = open(output_file, 'w') 

    print('\t'.join(["Chr", "Start", "End", "Ref", "Alt", "Original_Ref", "Variant_Ratio", "Depth", 
                     "Variant_Num", "Strand_Ratio", "Variant_Num_Info"]), file = hout)

    with open(input_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
        
            # first simple check
            total_variant_num = int(F[3]) - F[4].count('.') - F[4].count(',') - F[4].count('>') - F[4].count('<') + F[4].count('+') + F[4].count('-')
            if total_variant_num < min_variant_num: continue
            if total_variant_num / int(F[3]) < min_variant_ratio: continue

            var2num = {}
            var2pos = {}
            var2num_plus = {}
            bases = F[4]
            # qualities = F[5].split('')
            # positions = F[6].split(',')
            base_ind = 0
            depth_p, depth_n = 0, 0
            while bases != '':
                if bases[0] in ['>', '<', '*']: 
                    base_ind = base_ind + 1
                    bases = bases[1:]

                elif bases[0] in '^':
                    bases = bases[2:]
                elif bases[0] in '$':
                    bases = bases[1:]
                elif bases[0] in ['.', ',', 'A', 'C', 'G', 'T', 'N', 'a', 'c', 'g', 't', 'n']:
                    if bases[0] not in ['.', ',']: 
                        var_original = bases[0]
                        var = var_original.upper()
                        if var not in var2num:
                            var2num[var], var2pos[var], var2num_plus[var] = 0, [], 0
                        var2num[var] = var2num[var] + 1
                        # var2pos[var].append(positions[base_ind])
                        if var == var_original: 
                            var2num_plus[var] = var2num_plus[var] + 1

                    if bases[0] in ['.', 'A', 'C', 'G', 'T', 'N']:
                        depth_p = depth_p + 1
                    else:
                        depth_n = depth_n + 1

                    bases = bases[1:]

                    
                    if len(bases) > 0 and bases[0] in ['+', '-']:

                        match = re.search(r'^[\+\-](\d+)', bases)
                        indel_size = int(match.group(1))
                        # var_original = bases[0] + bases[2:(2 + indel_size)]
                        var_original = bases[0] + bases[(len(str(indel_size)) + 1):(len(str(indel_size)) + indel_size + 1)]
                        var = var_original.upper()
                        if var not in var2num:
                            var2num[var], var2pos[var], var2num_plus[var] = 0, [], 0
                        var2num[var] = var2num[var] + 1
                        # var2pos[var].append(positions[base_ind])
                        if var == var_original: var2num_plus[var] = var2num_plus[var] + 1

                        # bases = bases[(2 + indel_size):]
                        bases = bases[(len(str(indel_size)) + indel_size + 1):]
                    base_ind = base_ind + 1

            # if len(positions) != base_ind:
            #     print("Error???")
            #     sys.exit(1)

            if depth_p + depth_n == 0: continue


            bvar = ''
            bmis_rate = 0
            for var in var2num:
                if var2num[var] < min_variant_num: continue
                cur_rate = float(var2num[var]) / (depth_p + depth_n)
                if cur_rate > bmis_rate:
                    bmis_rate = cur_rate
                    bvar = var

            if bmis_rate < min_variant_ratio: continue
            
            # unique_positions = list(set([int(x) for x in var2pos[bvar]]))
            # if len(unique_positions) < min_variant_num: continue
            # if max(unique_positions) - min(unique_positions) < min_pos_range: continue
            

            var_info = str(depth_p) + ',' + str(var2num_plus[bvar]) + ';' + str(depth_n) + ',' + str(var2num[bvar] - var2num_plus[bvar])
            strand_ratio = float(var2num_plus[bvar]) / var2num[bvar]

            start, end, ref, alt, original_ref = F[1], F[1], F[2], bvar, F[2]
            if bvar.startswith('-'): 
                ref = bvar[1:]
                alt = '-'
                start, end = str(int(start) + 1), str(int(start) + len(ref))
            if bvar.startswith('+'):
                ref = '-'
                alt = bvar[1:]

            print('\t'.join([F[0], start, end, ref, alt, original_ref, str(round(bmis_rate, 4)), str(depth_p + depth_n), str(var2num[bvar]),
                             str(round(strand_ratio, 4)), var_info]), file = hout)


    hout.close()


def alu_mut_org(output_h, input_file, lid, lpos, lstrand):
    
    with open(input_file, 'r') as hin:
        header = hin.readline().rstrip('\n').split('\t')
        for line in hin:
            F = line.rstrip('\n').split('\t')

            if int(F[1]) <= 500 or int(F[2]) >= 812: continue            
            if lstrand == '+':
                F[1] = str(int(F[1]) - 500)
                F[2] = str(int(F[2]) - 500)
            else:
                temp_end = str(812 - int(F[1]))
                F[1] = str(812 - int(F[2]))
                F[2] = temp_end
                F[3] = reverse_complement(F[3]) if F[3] != "-" else "-"
                F[4] = reverse_complement(F[4]) if F[4] != "-" else "-"
                F[5] = reverse_complement(F[5]) if F[5] != "-" else "-"

            print('\t'.join([lid, lpos, lstrand] + F[1:]), file = output_h)



def alu_genotype_main(input_file, output_file, reference):

    def mut_call_proc(output_h, temp_key, tmp_dir, reference, reference_h):

        lchr, lpos, lid, lstrand, ltsd = temp_key.split(',')
        with open(tmp_dir + '/' + temp_key + ".ref.fa", 'w') as hout:
    
            print('>' + lid, file = hout)
            print(generate_alu_seq(reference_h, lchr, lpos, lstrand, ltsd), file = hout) 

        subprocess.check_call(["bwa", "index", tmp_dir + '/' + temp_key + ".ref.fa"], stderr = subprocess.DEVNULL)
        with open(tmp_dir + '/' + temp_key + ".sam", 'w') as hout:
            subprocess.check_call(["bwa", "mem", tmp_dir + '/' + temp_key + ".ref.fa", 
                                   tmp_dir + '/' + temp_key + "_1.fq", tmp_dir + '/' + temp_key + "_2.fq"], stdout = hout, stderr = subprocess.DEVNULL)

        with open(tmp_dir + '/' + temp_key + ".bam", 'w') as hout:
            subprocess.check_call(["samtools", "view", "-bhS", tmp_dir + '/' + temp_key + ".sam"], stdout = hout, stderr = subprocess.DEVNULL)
    
        with open(tmp_dir + '/' + temp_key + ".sorted.bam", 'w') as hout:
            subprocess.check_call(["samtools", "sort", tmp_dir + '/' + temp_key + ".bam"], stdout = hout, stderr = subprocess.DEVNULL)

        with open(tmp_dir + '/' + temp_key + ".pileup", 'w') as hout: 
            subprocess.check_call(["samtools", "mpileup", "-f", tmp_dir + '/' + temp_key + ".ref.fa", tmp_dir + '/' + temp_key + ".sorted.bam"], 
                                  stdout = hout, stderr = subprocess.DEVNULL)
 
        proc_mpileup(tmp_dir + '/' + temp_key + ".pileup", tmp_dir + '/' + temp_key + ".mut.txt")

        alu_mut_org(output_h, tmp_dir + '/' + temp_key + ".mut.txt", lid, lpos, lstrand)



    reference_h = pysam.FastaFile(os.path.abspath(reference))
    hout_m = open(output_file, 'w')

    tmp_dir = "tmp"
    p_tmp_dir = pathlib.Path(tmp_dir)
    if not p_tmp_dir.exists(): p_tmp_dir.mkdir()

    temp_key = '' 
    with open(input_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            # if F[0] == "chr1,23654652,ALU_umary_ALU_69,+,AACAGATGGGGCATC":
            if F[0] != temp_key:

                if temp_key != '':
                    hout1.close()
                    hout2.close()

                    mut_call_proc(hout_m, temp_key, tmp_dir, reference, reference_h)

                temp_key = F[0]
                hout1 = open(tmp_dir + '/' + temp_key + "_1.fq", 'w')
                hout2 = open(tmp_dir + '/' + temp_key + "_2.fq", 'w')

            if F[1].endswith('1'):
                print('@' + F[1] + '\n' + F[2] + '\n' + '+' + '\n' + F[3], file = hout1)
            else:
                print('@' + F[1] + '\n' + F[2] + '\n' + '+' + '\n' + F[3], file = hout2)


    if temp_key != 0:
        hout1.close()
        hout2.close()
        mut_call_proc(hout_m, temp_key, tmp_dir, reference, reference_h)

    reference_h.close()
    hout_m.close()
        


if __name__ == "__main__":

    import sys
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    reference = sys.argv[3]

    # reference_h = pysam.FastaFile(os.path.abspath(reference))

    alu_genotype_main(input_file, output_file, reference)
 
