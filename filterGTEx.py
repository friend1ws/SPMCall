#! /usr/local/bin/python

"""
    script for extracting junctions filtering those appearing in GTEx data
    junctions with more than $thres0 read is extrancted and 
    junctions appearing >= $thres1 reads in >= $thres2 GTEx samples are filtered
"""

import sys, tabix

inputFile = sys.argv[1]
GTExFile = sys.argv[2]
thres0 = sys.argv[3]
thres1 = sys.argv[4]
thres2 = sys.argv[5]

hIN = open(inputFile, 'r')
gtex_tb = tabix.open(GTExFile)

for line in hIN:
    F = line.rstrip('\n').split('\t')

    # check control data 
    tabixErrorFlag = 0
    try:
        records = gtex_tb.query(F[0], int(F[1]) - 2, int(F[1]) + 2)
    except Exception as inst:
        print >> sys.stderr, "%s: %s at the following key:" % (type(inst), inst.args)
        print >> sys.stderr, '\t'.join(F)
        tabixErrorFlag = 1

    gtexFlag = 0
    if tabixErrorFlag == 0:
        for record in records:
            if F[0] == record[0] and F[1] == record[1] and F[2] == record[2]:
                if len([x for x in record[3:] if int(x) >= int(thres1)]) >= int(thres2):
                    gtexFlag = 1

    if int(F[22]) >= int(thres0) and gtexFlag == 0:
        print '\t'.join(F[0:3] + F[22:23])

hIN.close()

