#! /usr/local/bin/python

import sys, tabix

inputFile = sys.argv[1]
controlFile = sys.argv[2]
thres = sys.argv[3]

hIN = open(inputFile, 'r')
control_tb = tabix.open(controlFile)

for line in hIN:
    F = line.rstrip('\n').split('\t')

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
                controlFlag = 1

    if int(F[22]) >= int(thres) and controlFlag == 0:
        print '\t'.join(F[0:3] + F[22:23])

hIN.close()

