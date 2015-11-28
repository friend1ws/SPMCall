#! /usr/bin/env python

import pysam
import datetime

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
            if not line.startswith('#'): print line

    hout.close()

