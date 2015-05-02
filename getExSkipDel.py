#! /usr/local/bin/python

import sys, tabix
from scipy import stats

inputFile = sys.argv[1]
tumorDepth = sys.argv[2]
normalDepth = sys.argv[3]

hIN = open(inputFile, 'r')
tumorDepth_tb = tabix.open(tumorDepth)
normalDepth_tb = tabix.open(normalDepth)

margin1 = 1000
margin2 = 500000
thres = 50

# for sorting
def cmp_chrPos(x1, x2):
    key1 = x1.split('\t')
    key2 = x2.split('\t')
    
    if key1[0] < key2[0]:
        return 1
    elif key1[0] > key2[0]:
        return -1
    else:
        if int(key1[1]) >= int(key2[1]):
            return 1
        else:
            return -1


for line in hIN:
    F = line.rstrip('\n').split('\t')
    
    if F[3] != "exon-skip": continue
    searchRegion1 = [F[0], int(F[1]) + 10, int(F[2]) - 10]
    searchRegion2 = [F[0], int(F[1]) - margin2, int(F[1]) - margin1]
    searchRegion3 = [F[0], int(F[2]) + margin1, int(F[2]) + margin2]

    targetDepth = {} 
    flankingDepth1 = {} 
    flankingDepth2 = {}

    # rough check for the depth between the spliced region
    tabixErrorFlag1 = 0
    try:
        tumorDepth1 = tumorDepth_tb.query(searchRegion1[0], searchRegion1[1], searchRegion1[2])
    except Exception as inst:
        print >> sys.stderr, "%s: %s at the following key:" % (type(inst), inst.args)
        print >> sys.stderr, '\t'.join(F)
        tabixErrorFlag1 = 1
 
    tabixErrorFlag2 = 0
    try:
        normalDepth1 = normalDepth_tb.query(searchRegion1[0], searchRegion1[1], searchRegion1[2])
    except Exception as inst:
        print >> sys.stderr, "%s: %s at the following key:" % (type(inst), inst.args)
        print >> sys.stderr, '\t'.join(F)
        tabixErrorFlag2 = 1

    if tabixErrorFlag1 == 0 and tabixErrorFlag2 == 0:
        tempKey2depth1_T = {}
        tempKey2depth1_N = {}
        for record in tumorDepth1:
            tempKey2depth1_T[record[0] + '\t' + record[1] + '\t' + record[2]] = int(record[3])
        for record in normalDepth1: 
            if int(record[3]) >= thres: tempKey2depth1_N[record[0] + '\t' + record[1] + '\t' + record[2]] = int(record[3])

        for key in tempKey2depth1_N:
            if key in tempKey2depth1_T: targetDepth[key] = float(tempKey2depth1_T[key]) / float(tempKey2depth1_N[key])


    if len(targetDepth) == 0: continue

    # rough check for the depth outside (left side) the spliced region
    tabixErrorFlag1 = 0
    try:
        tumorDepth2 = tumorDepth_tb.query(searchRegion2[0], searchRegion2[1], searchRegion2[2])
    except Exception as inst:
        print >> sys.stderr, "%s: %s at the following key:" % (type(inst), inst.args)
        print >> sys.stderr, '\t'.join(F)
        tabixErrorFlag1 = 1
    
    tabixErrorFlag2 = 0
    try:
        normalDepth2 = normalDepth_tb.query(searchRegion2[0], searchRegion2[1], searchRegion2[2])
    except Exception as inst:
        print >> sys.stderr, "%s: %s at the following key:" % (type(inst), inst.args)
        print >> sys.stderr, '\t'.join(F)
        tabixErrorFlag2 = 1
        
    if tabixErrorFlag1 == 0 and tabixErrorFlag2 == 0:
        tempKey2depth2_T = {}
        tempKey2depth2_N = {}
        for record in tumorDepth2:
            tempKey2depth2_T[record[0] + '\t' + record[1] + '\t' + record[2]] = int(record[3])
        for record in normalDepth2: 
            if int(record[3]) >= thres: tempKey2depth2_N[record[0] + '\t' + record[1] + '\t' + record[2]] = int(record[3])
            
        for key in tempKey2depth2_N:
            if key in tempKey2depth2_T: flankingDepth1[key] = float(tempKey2depth2_T[key]) / float(tempKey2depth2_N[key])


    # rough check for the depth outside (left side) the spliced region
    tabixErrorFlag1 = 0
    try:
        tumorDepth3 = tumorDepth_tb.query(searchRegion3[0], searchRegion3[1], searchRegion3[2])
    except Exception as inst:
        print >> sys.stderr, "%s: %s at the following key:" % (type(inst), inst.args)
        print >> sys.stderr, '\t'.join(F)
        tabixErrorFlag1 = 1

    tabixErrorFlag2 = 0
    try:
        normalDepth3 = normalDepth_tb.query(searchRegion3[0], searchRegion3[1], searchRegion3[2])
    except Exception as inst:
        print >> sys.stderr, "%s: %s at the following key:" % (type(inst), inst.args)
        print >> sys.stderr, '\t'.join(F)
        tabixErrorFlag2 = 1
 
    if tabixErrorFlag1 == 0 and tabixErrorFlag2 == 0:
        tempKey2depth2_T = {}
        tempKey2depth2_N = {}
        for record in tumorDepth3:
            tempKey2depth2_T[record[0] + '\t' + record[1] + '\t' + record[2]] = int(record[3])
        for record in normalDepth3:
            if int(record[3]) >= thres: tempKey2depth2_N[record[0] + '\t' + record[1] + '\t' + record[2]] = int(record[3])

        for key in tempKey2depth2_N:
            if key in tempKey2depth2_T: flankingDepth2[key] = float(tempKey2depth2_T[key]) / float(tempKey2depth2_N[key])

    

    flankingDepth = {}
    fLen1 = len(flankingDepth1)
    num = 1
    for key in sorted(flankingDepth1, cmp = cmp_chrPos): 
        if fLen1 - num <= 10:
            flankingDepth[key] = flankingDepth1[key]
        num = num + 1

    fLen2 = len(flankingDepth2)
    num = 1
    for key in sorted(flankingDepth2, cmp = cmp_chrPos):
        if fLen2 - num <= 10:
            flankingDepth[key] = flankingDepth2[key]
        num = num + 1


    wilcox_result = stats.ranksums(targetDepth.values(), flankingDepth.values())
    if float(wilcox_result[0]) > 0 or float(wilcox_result[1]) > 0.1: continue

    
    print "exon" + '\t' + '\t'.join(F[0:3]) + '\t' + str(wilcox_result[0]) + '\t' + str(wilcox_result[1])
    print "target"
    for key in sorted(targetDepth, cmp = cmp_chrPos):
        print key + '\t' + str(targetDepth[key]) + '\t' + str(tempKey2depth1_T[key]) + '\t' + str(tempKey2depth1_N[key])
    print "flanking"
    for key in sorted(flankingDepth, cmp = cmp_chrPos):
        # print key + '\t' + str(flankingDepth[key])
        print key + '\t' + str(flankingDepth[key]) + '\t' + str(tempKey2depth2_T[key]) + '\t' + str(tempKey2depth2_N[key])


hIN.close()
 

