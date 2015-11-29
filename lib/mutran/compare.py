#! /usr/bin/env python

import sys, tabix, numpy

def comp_mut_junction(input_file, output_file, mutation_file, annotation_dir):

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

    searchMargin1 = 30
    searchMargin2 = 10

    splicingDonnorMotif = ["AG", "GTRAGT"]
    splicingAcceptorMotif = ["YYYYNCAG", "G"]

    nuc2vec = {'A': [1, 0, 0, 0], 'C': [0, 1, 0, 0], 'G': [0, 0, 1, 0], 'T': [0, 0, 0, 1], \
               'W': [1, 0, 0, 1], 'S': [0, 1, 1, 0], 'M': [1, 1, 0, 0], 'K': [0, 0, 1, 1], \
               'R': [1, 0, 1, 0], 'Y': [0, 1, 0, 1], 'B': [0, 1, 1, 1], 'D': [1, 0, 1, 1], \
               'H': [1, 1, 0, 1], 'V': [1, 1, 1, 0], 'N': [1, 1, 1, 1]}

    ref_exon_bed = annotation_dir + "/refExon.bed.gz"
    grch2ucsc_file = annotation_dir + "/grch2ucsc.txt"

    # relationship between CRCh and UCSC chromosome names
    grch2ucsc = {}
    with open(grch2ucsc_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            grch2ucsc[F[0]] = F[1]

    hin = open(input_file, 'r')
    mutation_tb = tabix.open(mutation_file)
    exon_tb = tabix.open(ref_exon_bed)

    for line in hin:
        F = line.rstrip('\n').split('\t') 
        chr_name = grch2ucsc[F[0]] if F[0] in grch2ucsc else F[0]

        if F[3] not in ["exon-skip", "splice-site-slip", "pseudo-exon-inclusion"]: continue
        firstSearchRegion = [chr_name, int(F[1]), int(F[2])]
        splicingMotifRegions = []
        targetGene  =[]

        # we need to detect the non exon-intron junction break points
        # current procedure may be not perfect and be subject to change..
        gene1 = F[4].split(';')
        gene2 = F[7].split(';')
        junction1 = F[6].split(';')   
        junction2 = F[9].split(';')


        # just consider genes sharing the exon-intron junction with the breakpoints of splicings
        for i in range(0, len(gene1)):
            if junction1[i] != "*": targetGene.append(gene1[i])
        for i in range(0, len(gene2)):
            if junction2[i] != "*": targetGene.append(gene2[i])
        targetGene = list(set(targetGene))


        if F[3] in ["splice-site-slip", "pseudo-exon-inclusion"]:
            # for non exon-intron junction breakpoints
            if "*" in junction1 and "s" in junction2: # splicing donnor motif, plus direction
                firstSearchRegion[1] = firstSearchRegion[1] - searchMargin1
                splicingMotifRegions.append((chr_name, int(F[1]) - len(splicingDonnorMotif[0]) + 1, int(F[1]) + len(splicingDonnorMotif[1]), "donnor", "+", 0))
            if "*" in junction1 and "e" in junction2: # splicing acceptor motif, minus direction
                firstSearchRegion[1] = firstSearchRegion[1] - searchMargin1
                splicingMotifRegions.append((chr_name, int(F[1]) - len(splicingAcceptorMotif[1]) + 1, int(F[1]) + len(splicingAcceptorMotif[0]), "acceptor", "-", 0))
            if "s" in junction1 and "*" in junction2: # splicing donnor motif, minus direction
                firstSearchRegion[2] = firstSearchRegion[2] + searchMargin1
                splicingMotifRegions.append((chr_name, int(F[2]) - len(splicingDonnorMotif[1]), int(F[2]) + len(splicingDonnorMotif[0]) - 1, "donnor", "-", 0))
            if "e" in junction1 and "*" in junction2: # # splicing acceptor motif, plus direction
                firstSearchRegion[2] = firstSearchRegion[2] + searchMargin1
                splicingMotifRegions.append((chr_name, int(F[2]) - len(splicingAcceptorMotif[0]), int(F[2]) + len(splicingAcceptorMotif[1]) - 1, "acceptor", "+", 0))


        ##########
        # rough check for the mutation between the spliced region
        tabixErrorFlag1 = 0
        try:
            mutations = mutation_tb.query(firstSearchRegion[0], firstSearchRegion[1], firstSearchRegion[2])
        except Exception as inst:
            # print >> sys.stderr, "%s: %s at the following key:" % (type(inst), inst.args)
            # print >> sys.stderr, '\t'.join(F)
            tabixErrorFlag1 = 1

        # if there are some mutaions
        if tabixErrorFlag1 == 0 and mutations is not None:

            # check the exons within the spliced regions
            tabixErrorFlag2 = 0
            try:
                exons = exon_tb.query(firstSearchRegion[0], firstSearchRegion[1], firstSearchRegion[2])
            except Exception as inst:
                # print >> sys.stderr, "%s: %s at the following key:" % (type(inst), inst.args)
                # print >> sys.stderr, '\t'.join(F)
                tabixErrorFlag2 = 1

            # first, add the exon-intron junction for detailed check region list
            if tabixErrorFlag2 == 0:
                for exon in exons:
                    if exon[3] not in targetGene: continue
                    if exon[5] == "+":
                        # splicing acceptor motif, plus direction
                        splicingMotifRegions.append((exon[0], int(exon[1]) - len(splicingAcceptorMotif[0]) + 1, int(exon[1]) + len(splicingAcceptorMotif[1]), "acceptor", "+", 1))
                        # splicing donnor motif, plus direction
                        splicingMotifRegions.append((exon[0], int(exon[2]) - len(splicingDonnorMotif[0]) + 1, int(exon[2]) + len(splicingDonnorMotif[1]), "donnor", "+", 1))
                    if exon[5] == "-":
                        # splicing donnor motif, minus direction 
                        splicingMotifRegions.append((exon[0], int(exon[1]) - len(splicingDonnorMotif[1]) + 1, int(exon[1]) + len(splicingDonnorMotif[0]), "donnor", "-", 1))
                        # splicing acceptor motif, minus direction
                        splicingMotifRegions.append((exon[0], int(exon[2]) - len(splicingAcceptorMotif[1]) + 1, int(exon[2]) + len(splicingAcceptorMotif[0]), "acceptor", "-", 1))


            splicingMotifRegions = list(set(splicingMotifRegions))

            if F[1] == "56489094" and F[2] == "56490287":
                pass

            # compare each mutation with exon-intron junction regions and non-exon-intorn junction breakpoints.
            for mutation in mutations:
                RegMut = []
                for reg in splicingMotifRegions:

                    # insertion or deletion (just consider the disruption of splicing motifs now)
                    if (len(mutation[3]) > 1 or len(mutation[4]) > 1) and reg[5] == 1:
                        if int(mutation[1]) <= reg[2] - 1 and reg[1] <= int(mutation[1]) + len(mutation[4]) - 1:
                            RegMut.append([reg, "splicing " + reg[3] + " disruption"])
                    
                    # base substitution
                    if len(mutation[3]) == 1 and len(mutation[4]) == 1 and reg[1] <= int(mutation[1]) <= reg[2]:
                        motifSeq = ""
                        if reg[3] == "acceptor": motifSeq = splicingAcceptorMotif[0] + splicingAcceptorMotif[1]                  
                        if reg[3] == "donnor": motifSeq = splicingDonnorMotif[0] + splicingDonnorMotif[1]

                        if reg[4] == "-":
                            complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', \
                                          'W': 'S', 'S': 'W', 'M': 'K', 'K': 'M', \
                                          'R': 'Y', 'Y': 'R', 'B': 'V', 'D': 'H', \
                                          'H': 'D', 'V': 'B', 'N': 'N'}
                            motifSeq = "".join(complement.get(base) for base in reversed(motifSeq)) 


                        vecAtMut = nuc2vec[motifSeq[int(mutation[1]) - int(reg[1])]]
                        editDistDiff = numpy.dot(vecAtMut, nuc2vec[mutation[4]]) - numpy.dot(vecAtMut, nuc2vec[mutation[3]]) 

                        if editDistDiff > 0 and reg[5] == 0: RegMut.append([reg, "splicing " + reg[3] + " creation"])
                        if editDistDiff < 0 and reg[5] == 1: RegMut.append([reg, "splicing " + reg[3] + " disruption"])

                for item in RegMut:
                    print >> hout, '\t'.join(F) + '\t' + '\t'.join(mutation) + '\t' + item[0][0] + ':' + str(item[0][1]) + '-' + str(item[0][2]) + ',' + item[0][4] + '\t' + item[1]


    hin.close()
    hout.close()

