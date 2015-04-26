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
