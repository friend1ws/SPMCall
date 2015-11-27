#! /usr/local/bin/python

import sys, tabix

inputFile = sys.argv[1]
targetFile = sys.argv[2]
thres = sys.argv[3]

target_tb = tabix.open(targetFile)
hIN = open(inputFile, 'r')
for line in hIN:
    
    F = line.rstrip('\n').split('\t')
    if int(F[22]) < int(thres): continue

    # check target data 
    tabixErrorFlag = 0
    try:
        records = target_tb.query(F[0], int(F[1]) - 2, int(F[1]) + 2)
    except Exception as inst:
       print >> sys.stderr, "%s: %s at the following key:" % (type(inst), inst.args)
       print >> sys.stderr, '\t'.join(F)
       tabixErrorFlag = 1

    targetFlag = 0
    if tabixErrorFlag == 0:
        for record in records:
            if F[0] == record[0] and F[1] == record[1] and F[2] == record[2]:
                targetFlag = 1 

    
    if targetFlag == 1:
        print '\t'.join(F[0:3] + F[22:23])

hIN.close()
 
