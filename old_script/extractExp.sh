#! /bin/sh
#$ -S /bin/sh

TISSUE=$1

echo "python extractExp.py../data/GTEx/GTEx_Analysis_2014-01-17_RNA-seq_Flux1.6_junction_reads.txt.gz ../data/GTEx/GTEx_Data_2014-01-17_Annotations_SampleAttributesDS.txt ${TISSUE} > ../data/GTEx/${TISSUE}.junctions.txt"
python extractExp.py ../data/GTEx/GTEx_Analysis_2014-01-17_RNA-seq_Flux1.6_junction_reads.txt.gz ../data/GTEx/GTEx_Data_2014-01-17_Annotations_SampleAttributesDS.txt ${TISSUE} > ../data/GTEx/${TISSUE}.junctions.txt

echo "sort -k1,1 -k2,2n -k3,3 ../data/GTEx/${TISSUE}.junctions.txt > ../data/GTEx/${TISSUE}.junctions.bed"
sort -k1,1 -k2,2n -k3,3 ../data/GTEx/${TISSUE}.junctions.txt > ../data/GTEx/${TISSUE}.junctions.bed

echo "bgzip -f ../data/GTEx/${TISSUE}.junctions.bed > ../data/GTEx/${TISSUE}.junctions.bed.gz"
bgzip -f ../data/GTEx/${TISSUE}.junctions.bed > ../data/GTEx/${TISSUE}.junctions.bed.gz

echo "tabix -p bed -f ../data/GTEx/${TISSUE}.junctions.bed.gz"
tabix -p bed -f ../data/GTEx/${TISSUE}.junctions.bed.gz

rm -rf ../data/GTEx/${TISSUE}.junctions.txt


