#! /usr/local/bin/python

"""
    The purpose of this script is to classify splicing changes
    mainly by comparing the two breakpoints with the exon-intorn junction of genes
    within the database.
    Also, we generate the sequence arrond the breakpoints, which will be helpful
    for checking the authenticity of the splicing and evaluating the relationships
    with the somatic mutations.

    here is the classification categories:
    1. known (The splicing pattern is included in the database)
    (the start and end breakpoints are next exon-intron junctions of the same gene) 
    2. exon skipping 
    (the start and end breakpoints are exon-intron junctions of the same gene,
     but not the next ones)
    3. splice-site slip
    (one of the two breakpoints is an exon-intron junction and the other is within the 30bp exon of the same gene)
    4. pseudo-exon inclusion
    (one of the two break points is an exon-intron junction and the other is located in the same gene, but more than 30bp from exons of the gene)
    5. other
    (neighter of the two breakpoins are exon-intron junction, but located in the same gene)
    6. chimeric (spliced)
    7. chimeric (un-spliced)


    The algorithm for the annotation is as follows
    1. for both breakpoints, list up the exon-intron junctions matching to the breakpoints
    2. for both breakpoints, list up the exons within 30bp from the breakpoints
    3. for both breakpoints, list up the genes matching to the breakpoints
    4. summarize the above results and induce the annotation from them
    5. get the sequence arround the breakpoints.

"""

import sys, tabix

inputFile = sys.argv[1]
geneFile = sys.argv[2]
exonFile = sys.argv[3]

junction_margin = 3
exon_margin = 30

hIN = open(inputFile, 'r')
gene_tb = tabix.open(geneFile)
exon_tb = tabix.open(exonFile)

for line in hIN:
    F = line.rstrip('\n').split('\t')

    ##########
    # check gene annotation for the side 1  
    tabixErrorFlag = 0
    try:
        records = gene_tb.query(F[0], int(F[1]) - 1, int(F[1]))
    except Exception as inst:
        print >> sys.stderr, "%s: %s at the following key:" % (type(inst), inst.args)
        print >> sys.stderr, '\t'.join(F)
        tabixErrorFlag = 1

    gene1 = [];
    if tabixErrorFlag == 0:
        for record in records:
            gene1.append(record[3])

    gene1 = list(set(gene1))
    ##########

    ##########
    # check gene annotation for the side 2  
    tabixErrorFlag = 0
    try:
        records = gene_tb.query(F[0], int(F[2]) - 1, int(F[2]))
    except Exception as inst:
        print >> sys.stderr, "%s: %s at the following key:" % (type(inst), inst.args)
        print >> sys.stderr, '\t'.join(F)
        tabixErrorFlag = 1
        
    gene2 = [];
    if tabixErrorFlag == 0:
        for record in records:
            gene2.append(record[3])
            
    gene2 = list(set(gene2))
    ##########

    ##########
    # check exon and junction annotation for the side 1  
    tabixErrorFlag = 0
    try:
        records = exon_tb.query(F[0], int(F[1]) - exon_margin, int(F[1]) + exon_margin)
    except Exception as inst:
        print >> sys.stderr, "%s: %s at the following key:" % (type(inst), inst.args)
        print >> sys.stderr, '\t'.join(F)
        tabixErrorFlag = 1
    
    exon1 = {};
    junction1 = {};
    if tabixErrorFlag == 0:
        for record in records:
            exon1[record[3]] = record[4]
            if abs(int(F[1]) - int(record[2])) < junction_margin:
                if record[5] == "+": junction1[record[3]] = "end"
                if record[5] == "-": junction1[record[3]] = "start"
    ##########

    ##########
    # check exon and junction annotation for the side 2
    tabixErrorFlag = 0
    try:
        records = exon_tb.query(F[0], int(F[2]) - exon_margin, int(F[2]) + exon_margin)
    except Exception as inst:
        print >> sys.stderr, "%s: %s at the following key:" % (type(inst), inst.args)
        print >> sys.stderr, '\t'.join(F)
        tabixErrorFlag = 1

    exon2 = {};
    junction2 = {};
    if tabixErrorFlag == 0:
        for record in records:
            exon2[record[3]] = record[4]
            if abs(int(F[2]) - int(record[1])) < junction_margin:
                if record[5] == "+": junction2[record[3]] = "start"
                if record[5] == "-": junction2[record[3]] = "end"
    ##########


    spliceClass = ""
    checkGenes = list(set([gene1, gene2]))
    ##########
    # check for know junction
    passGene = []
    for gene in checkGenes:
        if gene in gene1 and gene in gene2 and gene in junction1 and gene in junction2:
            if junction1[gene] == "end" and junction2[gene] == "start" and exon2[gene] - exon1[gene] == 1: passGene.append(gene)
            if junction2[gene] == "end" and junction1[gene] == "start" and exon1[gene] - exon2[gene] == 1: passGene.append(gene)

    if len(passGene) > 1: spliceClass = "known"


    ##########
    # check for exon skip
    if spliceClass != "":
        passGene = []
        for gene in checkGenes:
            if gene in gene1 and gene in gene2 and gene in junction1 and gene in junction2:
                if junction1[gene] == "end" and junction2[gene] == "start" and exon2[gene] - exon1[gene] > 1: passGene.append(gene)
                if junction2[gene] == "end" and junction1[gene] == "start" and exon1[gene] - exon2[gene] > 1: passGene.append(gene)

    if len(passGene) > 1: spliceClass = "exon-skip"


    ##########
    # check for splice-site slip 
    if spliceClass != "":
        passGene = []
        for gene in checkGenes:
            if gene in gene1 and gene in gene2 and gene:
                if gene in junction1 and gene in exon2 and gene not in junction2: passGene.append(gene)
                if gene in junction2 and gene in exon1 and gene not in junction1: passGene.append(gene)

    if len(passGene) > 1: spliceClass = "splice-site-slip"
 

    ##########
    # check for pseudo-exon inclusion 
    if spliceClass != "":
        passGene = []
        for gene in checkGenes:
            if gene in gene1 and gene in gene2 and gene: 
                if gene in junction1 and gene not in exon2: passGene.append(gene)
                if gene in junction2 and gene not in exon1: passGene.append(gene)
                 
    if len(passGene) > 1: spliceClass = "pseudo-exon inclusion"


    ##########
    # check for within-gene 
    if spliceClass != "":
        passGene = []
        for gene in checkGenes:
            if gene in gene1 and gene in gene2 and gene: passGene.append(gene)   
    
    if len(passGene) > 1: spliceClass = "within-gene"


    ##########
    # check for spliced-chimera 
    if spliceClass != "":
        passGene = []
        for g1 in gene1:
            for g2 in gene2:
                if junction1[g1] == "start" and junction2[g2] == "end": passGene.append(g1 + ',' + g2)
                if junction1[g1] == "end" and junction2[g2] == "start": passGene.append(g1 + ',' + g2)

    if len(passGene) > 1: spliceClass = "spliced-chimera"


    ##########
    # check for unspliced-chimera 
    if spliceClass != "":
        passGene = []
        for g1 in gene1:
            for g2 in gene2:
                passGene.append(g1 + ',' + g2)

    if len(passGene) > 1: spliceClass = "unspliced-chimera"







    print '\t'.join(F[0:3]) + '\t' + '\t'.join([';'.join(gene1), ';'.join(exon1), ';'.join(junction1), ';'.join(gene2), ';'.join(exon2), ';'.join(junction2)]) 


hIN.close() 
