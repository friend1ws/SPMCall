#! /usr/local/bin/python

"""
    Script for merging the junction files in the "fileList", removing those registered in the GTEx data.
    Here the filter condition is 
    1. At least, one of the junction has more than "thres0" reads.
    2. Junctions appearing >= $thres1 reads in >= $thres2 GTEx samples are filtered.
"""

import sys, tabix

fileList = sys.argv[1]
controlFile = sys.argv[2]
thres0 = sys.argv[3]
thres1 = sys.argv[4]
thres2 = sys.argv[5]

sample2file = {}
hIN = open(fileList, 'r')
for line in hIN:
    F = line.rstrip('\n').split('\t')
    sample2file[F[0]] = F[1]
hIN.close()


junc2exp = {}
sampleNum = len(sample2file)
   
sampleInd = 0 
for sample in sorted(sample2file):

    print >> sys.stderr, sample
    hIN = open(sample2file[sample], 'r')
    for line in hIN:
        F = line.rstrip('\n').split('\t')

        junc = '\t'.join(F[0:3])

        if junc not in junc2exp:
            junc2exp[junc] = [0] * sampleNum 

        junc2exp[junc][sampleInd] = F[22]

    hIN.close()
    sampleInd = sampleInd + 1

control_tb = tabix.open(controlFile)
for junc in sorted(junc2exp):

    if max(junc2exp[junc]) >= int(thres0):

        F = junc.split('\t')
        # check control data 
        tabixErrorFlag = 0
        try:
            records = control_tb.query(F[0], int(F[1]) - 1, int(F[2]) + 1)
        except Exception as inst:
            print >> sys.stderr, "%s: %s at the following key:" % (type(inst), inst.args)
            print >> sys.stderr, '\t'.join(F)
            tabixErrorFlag = 1

        controlFlag = 0
        if tabixErrorFlag == 0:
            for record in records:
                if F[0] == record[0] and F[1] == record[1] and F[2] == record[2]:
                    if len([x for x in record[3:] if int(x) >= int(thres1)]) >= int(thres2):
                        controlFlag = 1

        if controlFlag == 0:
            print junc + '\t' + '\t'.join([str(junc2exp[junc][i]) for i in range(sampleNum)])
 
