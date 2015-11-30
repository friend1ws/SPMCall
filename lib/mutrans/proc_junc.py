#! /usr/bin/env python

def proc_star_junction(input_file, output_file, read_num_thres, overhang_thres, remove_annotated):
    
    if read_num_thres is None: read_num_thres = 0
    if overhang_thres is None: overhang_thres = 0
    if remove_annotated is None: remove_annotated = False
    
    hout = open(output_file, 'w')
    with open(input_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            if remove_annotated == True and F[5] == "0": continue
            if int(F[6]) < read_num_thres: continue
            if int(F[8]) < overhang_thres: continue

            print >> hout, '\t'.join(F)

    hout.close()


