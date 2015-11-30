#! /usr/bin/env python

import subprocess


def make_control_junction(input_file_list, output_file, read_num_thres, overhang_thres, remove_annotated):

    """
    script for generating control junction database from the list of star junction (.SJ.tab.out) files
    """

    hin = open(input_file_list, 'r')
    hout = open(output_file + ".unsroted", 'w')

    with open(input_file_list, 'r') as hin:
        for line in hin:

            junction_file = line.rstrip('\n')
            with open(junction_file, 'r') as hin2:
                for line2 in hin2:

                    F = line2.rstrip('\n').split('\t')
                    if remove_annotated == True and F[5] != "0": continue
                    if int(F[6]) < read_num_thres: continue 
                    if int(F[8]) < overhang_thres: continue 
         
                    # convert to map-splice2 coordinate
                    # F[1] = str(int(F[1]) - 1)
                    # F[2] = str(int(F[2]) + 1)
         
                    print >> hout, '\t'.join(F)
                

    hout = open(output_file + ".sorted", 'w')
    s_ret = subprocess.call(["sort", "-k1,1", "-k2,2n", "-k3,3n", output_file + ".unsroted"], stdout = hout)
    hout.close()

    if s_ret != 0:
        print >> sys.stderr, "Error in sorting merged junction file"
        sys.exit(1)


    hout = open(output_file + ".merged", 'w')
    with open(output_file + ".sorted", 'r') as hin:
        temp_key = ""
        temp_read_num = 0
        for line in hin:
            F = line.rstrip('\n').split('\t')
            key = F[0] + '\t' + F[1] + '\t' + F[2]
            read_num = int(F[6])
            if key != temp_key:
                if temp_key != "":
                    print >> hout, temp_key + '\t' + str(read_num)
                temp_key = key
            else:
                if read_num > temp_read_num:
                    temp_read_num = read_num

        print >> hout, temp_key + '\t' + str(read_num)


    hout = open(output_file, 'w')
    s_ret = subprocess.call(["bgzip", "-f", "-c", output_file + ".merged"], stdout = hout)
    hout.close()

    if s_ret != 0:
        print >> sys.stderr, "Error in compression merged junction file"
        sys.exit(1)


    s_ret = subprocess.call(["tabix", "-p", "bed", output_file])
    if s_ret != 0:
        print >> sys.stderr, "Error in indexing merged junction file"
        sys.exit(1)


if __name__ == "__main__":

    import sys

    input_file_list = sys.argv[1]
    output_file = sys.argv[2]
    read_num_thres = int(sys.argv[3])
    overhang_thres = int(sys.argv[4])
    remove_annotated = bool(sys.argv[5])

    make_control_junction(input_file_list, output_file, read_num_thres, overhang_thres, remove_annotated)    
    

