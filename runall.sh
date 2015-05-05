#! /bin/sh
#$ -S /bin/sh
#$ -cwd

if [ ! -d ../result/ATL/annot ];
then
    mkdir -p ../result/ATL/filt
    mkdir -p ../result/ATL/annot
    mkdir -p ../result/ATL/mut2sp_WGS
fi

# filtering using GTEx data and negative binomial outlier test

# echo "python mergeAndFiltGTEx.py ../data/ATL/sample2juncPath.txt ../data/GTEx/Blood.junctions.bed.gz 5 3 2 > ../result/ATL/filt/merge.GTEx.filt.txt"
# python mergeAndFiltGTEx.py ../data/ATL/sample2juncPath.txt ../data/GTEx/Blood.junctions.bed.gz 5 3 2 > ../result/ATL/filt/merge.GTEx.filt.txt

# echo "Rscript nbDistFilt.R ../result/ATL/filt/merge.GTEx.filt.txt ../result/ATL/filt/merge.GTEx.filt.nb.txt"
# Rscript nbDistFilt.R ../result/ATL/filt/merge.GTEx.filt.txt ../result/ATL/filt/merge.GTEx.filt.nb.txt

# echo "python annotSplicing.py ../result/ATL/filt/merge.GTEx.filt.nb.txt db/refGene.bed.gz db/refExon.bed.gz | grep "exon-skip\|pseudo-exon-inclusion\|splice-site-slip" > ../result/ATL/filt/merge.GTEx.filt.nb.annot.txt"
# python annotSplicing.py ../result/ATL/filt/merge.GTEx.filt.nb.txt db/refGene.bed.gz db/refExon.bed.gz | grep "exon-skip\|pseudo-exon-inclusion\|splice-site-slip" > ../result/ATL/filt/merge.GTEx.filt.nb.annot.txt

# echo "sort -k1,1 -k2,2n -k3,3n ../result/ATL/filt/merge.GTEx.filt.nb.annot.txt | cut -f 1-3 | bgzip -f > ../result/ATL/filt/merge.GTEx.filt.nb.annot.bed.gz"
# sort -k1,1 -k2,2n -k3,3n ../result/ATL/filt/merge.GTEx.filt.nb.annot.txt | cut -f 1-3 | bgzip -f > ../result/ATL/filt/merge.GTEx.filt.nb.annot.bed.gz

# echo "tabix -p bed ../result/ATL/filt/merge.GTEx.filt.nb.annot.bed.gz"
# tabix -p bed ../result/ATL/filt/merge.GTEx.filt.nb.annot.bed.gz
 
# :<<_COMMENT_OUT_

while read sample;
do
    echo "python checkGTExNB.py ../data/ATL/junction/${sample}_T.junctions.txt ../result/ATL/filt/merge.GTEx.filt.nb.annot.bed.gz 5 > ../result/ATL/filt/${sample}.junctions.filt.txt"
    python checkGTExNB.py ../data/ATL/junction/${sample}_T.junctions.txt ../result/ATL/filt/merge.GTEx.filt.nb.annot.bed.gz 5 > ../result/ATL/filt/${sample}.junctions.filt.txt

    echo "python annotSplicing.py ../result/ATL/filt/${sample}.junctions.filt.txt db/refGene.bed.gz db/refExon.bed.gz > ../result/ATL/annot/${sample}.junctions.annot.txt"
    python annotSplicing.py ../result/ATL/filt/${sample}.junctions.filt.txt db/refGene.bed.gz db/refExon.bed.gz > ../result/ATL/annot/${sample}.junctions.annot.txt

done < ../data/ATL/sampleList_depth.txt
 
:<<_COMMENT_OUT_

while read sample1;
do
    while read sample2;
    do
        echo "python getSplicingMut.py ../result/ATL/annot/${sample1}.junctions.annot.txt ../data/ATL/mutation_WGS/${sample2}.mutation.vcf.gz db/refExon.bed.gz > ../result/ATL/mut2sp_WGS/${sample1}-${sample2}.mutation2splicing.txt"
        python getSplicingMut.py ../result/ATL/annot/${sample1}.junctions.annot.txt ../data/ATL/mutation_WGS/${sample2}.mutation.vcf.gz db/refExon.bed.gz > ../result/ATL/mut2sp_WGS/${sample1}-${sample2}.mutation2splicing.txt
    done < ../data/ATL/sampleList_WGS.txt
done < ../data/ATL/sampleList_WGS.txt

_COMMENT_OUT_

