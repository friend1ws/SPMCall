#! /usr/local/bin/python

"""
    a script for detecting candidate somatic substitutions causing splicing changes
    
    1. exon-skip:
       exon-intron junction region within spliced region
    
    2. splice-site-slip, pseudo-exon-inclusion
       in addition to exon-intron junction region, 
       we check for the region around non exon-intron junction break points 

    Maybe, we should check the exon-intron junction region within, 
    e.g., +-10000bp from the spliced region to search for unknown phenomena.

"""

import sys, tabix

spliceFile = sys.argv[1]
mutationFile = sys.argv[2]
exonFile = sys.argv[3]

searchMargin1 = 30
searchMargin2 = 10

hIN = open(spliceFile, 'r')
mutation_tb = tabix.open(mutationFile)
exon_tb = tabix.open(exonFile)

for line in hIN:
    F = line.rstrip('\n').split('\t') 
    if F[3] not in ["exon-skip", "splice-site-slip", "pseudo-exon-inclusion"]: continue
    firstSearchRegion = [F[0], int(F[1]), int(F[2])]
    detaiedSearchRegion = []

    # we need to detect the non exon-intron junction break points
    # current procedure may be not perfect and be subject to change..
    junction1 = F[6].split(';')   
    junction2 = F[9].split(';')

    if F[3] in ["splice-site-slip", "pseudo-exon-inclusion"]:
        if "*" in junction1:
            firstSearchRegion[1] = firstSearchRegion[1] - searchMargin1
            detaiedSearchRegion.append((F[0], int(F[1]) - searchMargin2, int(F[1]) + searchMargin2))
        if "*" in junction2:
            firstSearchRegion[2] = firstSearchRegion[2] + searchMargin1
            detaiedSearchRegion.append((F[0], int(F[2]) - searchMargin2, int(F[2]) + searchMargin2))


    ##########
    # rough check for the mutation between the spliced region
    tabixErrorFlag1 = 0
    try:
        mutations = mutation_tb.query(firstSearchRegion[0], firstSearchRegion[1], firstSearchRegion[2])
    except Exception as inst:
        print >> sys.stderr, "%s: %s at the following key:" % (type(inst), inst.args)
        print >> sys.stderr, '\t'.join(F)
        tabixErrorFlag1 = 1

    # if there are some mutaions
    if tabixErrorFlag1 == 0 and mutations is not None:

        # check the exons within the spliced regions
        tabixErrorFlag2 = 0
        try:
            exons = exon_tb.query(firstSearchRegion[0], firstSearchRegion[1], firstSearchRegion[2])
        except Exception as inst:
            print >> sys.stderr, "%s: %s at the following key:" % (type(inst), inst.args)
            print >> sys.stderr, '\t'.join(F)
            tabixErrorFlag2 = 1

        # first, add the exon-intron junction for detailed check region list
        if tabixErrorFlag2 == 0:
            for exon in exons:
                detaiedSearchRegion.append((exon[0], int(exon[1]) - searchMargin2, int(exon[1]) + searchMargin2))
                detaiedSearchRegion.append((exon[0], int(exon[2]) - searchMargin2, int(exon[2]) + searchMargin2))

        detaiedSearchRegion = list(set(detaiedSearchRegion))

        # compare each mutation with exon-intron junction regions and non-exon-intorn junction breakpoints.
        for mutation in mutations:
            for reg in detaiedSearchRegion:
                if reg[1] <= int(mutation[1]) <= reg[2]:
                    print '\t'.join(F) + '\t' + '\t'.join(mutation)

hIN.close()

