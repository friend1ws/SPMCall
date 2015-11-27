#! /usr/local/bin/python

"""
    script for extracting the junctions for specified tissue
"""

import sys, gzip

junctionFile = sys.argv[1] # "GTEx_Analysis_2014-01-17_RNA-seq_Flux1.6_junction_reads.txt.gz"
annotationFile = sys.argv[2] # "GTEx_Data_2014-01-17_Annotations_SampleAttributesDS.txt"
tissue = sys.argv[3]

targetIDs = []
hIN = open(annotationFile, 'r')
for line in hIN:
    F = line.rstrip('\n').split('\t')
    if F[0] == "SAMPID": continue # skip header

    if F[5] == tissue:
        targetIDs.append(F[0])

hIN.close()

 
targetInd = []
hIN = gzip.open(junctionFile, 'r')
for line in hIN:
    F = line.rstrip('\n').split('\t')

    if F[0] == "TargetID":
        for i in range(4, len(F)):
            if F[i] in targetIDs:
                targetInd.append(i)
    else:
        print "chr" + '\t'.join(F[0].split("_")) + '\t' + '\t'.join([str(int(float(F[i]))) for i in targetInd])

hIN.close()
 

    
