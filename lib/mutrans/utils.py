#! /usr/bin/env python

import pysam
import datetime
import tabix

def convert_anno2vcf(input_file, output_file, reference):

    """
        convert annovar format to vcf format
        the difficulty is: for the annovar format, the reference base for insertion or deletion is removed,
        and we have to recover them...
    """

    hin = open(input_file, 'r')
    hout = open(output_file, 'w')

    """
    # get current date
    today = datetime.date.today()

    # print meta-information lines
    print >> hout, "##fileformat=#VFCv4.1"
    print >> hout, "##fileDate=" + today.strftime("%Y%m%d")
    # print "##source=ebcallv2.0"
    # print "##INFO=<ID=TD,Number=1,Type=Integer,Description=\"Tumor Depth\">"
    # print "##INFO=<ID=TV,Number=1,Type=Integer,Description=\"Tumor Variant Read Num\">"
    # print "##INFO=<ID=ND,Number=1,Type=Integer,Description=\"Normal Depth\">"
    # print "##INFO=<ID=NV,Number=1,Type=Integer,Description=\"Normal Variant Read Num\">"
    print >> hout, "##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description=\"Somatic Variation\">"
    print >> hout, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
    """

    for line in hin:
        F = line.rstrip('\n').split('\t')
        
        pos, ref, alt = F[1], F[3], F[4]
        
        # insertion
        if F[3] == "-":
            # get the sequence for the reference base
            seq = ""    
            for item in pysam.faidx(reference, F[0] + ":" + str(F[1]) + "-" + str(F[1])):
                if item[0] == ">": continue
                seq = seq + item.rstrip('\n')
            ref, alt = seq, seq + F[4]

        # deletion
        if F[4] == "-":
            # get the sequence for the reference base
            seq = ""    
            for item in pysam.faidx(reference, F[0] + ":" + str(int(F[1]) - 1) + "-" + str(int(F[1]) - 1)):
                if item[0] == ">": continue
                seq = seq + item.rstrip('\n')

            pos, ref, alt = str(int(F[1]) - 1), seq + F[3], seq


        # QUAL = int(float(F[15]) * 10)
        QUAL = 60
        # INFO = "TD=" + F[9] + ";TV=" + F[10] + ";ND=" + F[13] + ";NV=" + F[14] + ";SOMATIC"
        INFO = "SOMATIC"

        print >> hout, F[0] + "\t" + pos + "\t.\t" + ref + "\t" + alt \
            + "\t" + str(QUAL) + "\t" + "PASS" + "\t" + INFO 



     
def remove_vcf_header(input_file, output_file):

    hout = open(output_file, 'w')
    with open(input_file, 'r') as hin:
        for line in hin:
            line = line.rstrip('\n')
            if not line.startswith('#'): print >> hout, line

    hout.close()


def proc_star_junction(input_file, output_file, control_file, read_num_thres, overhang_thres, remove_annotated):
    
    is_control = True if control_file is not None else False
    if is_control: control_db = tabix.open(control_file)

    if read_num_thres is None: read_num_thres = 0
    if overhang_thres is None: overhang_thres = 0
    if remove_annotated is None: remove_annotated = False
    
    hout = open(output_file, 'w')
    with open(input_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            key = F[0] + '\t' + F[1] + '\t' + F[2]
            if remove_annotated == True and F[5] != "0": continue
            if int(F[6]) < read_num_thres: continue
            if int(F[8]) < overhang_thres: continue

            if F[1].startswith("2959542"):
                pass

            ##########
            # remove control files
            if is_control:
                control_flag = 0
                tabixErrorFlag = 0
                try:
                    records = control_db.query(F[0], int(F[1]) - 5, int(F[1]) + 5)
                except Exception as inst:
                    tabixErrorFlag = 1

                control_flag = 0;
                if tabixErrorFlag == 0:
                    for record in records:
                        if F[0] == record[0] and F[1] == record[1] and F[2] == record[2]:
                            control_flag = 1

                if control_flag == 1: continue
            ##########

            # convert to map-splice2 coordinate
            F[1] = str(int(F[1]) - 1)
            F[2] = str(int(F[2]) + 1)

            print >> hout, '\t'.join(F)

    hout.close()


def convert_genosv2bed(input_file, output_file):

    hout = open(output_file, 'w')
    num = 1
    with open(input_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            chr1, chr2 = F[0], F[3]
            start1, end1 = str(int(F[1]) - 1), F[1]
            start2, end2 = str(int(F[4]) - 1), F[4]
            dir1, dir2 = F[2], F[5]
            name = "SV_" + str(num)
            score = "0"

            print >> hout, '\t'.join([chr1, start1, end1, chr2, start2, end2, name, score, dir1, dir2])

    hout.close()

