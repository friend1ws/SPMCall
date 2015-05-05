#! /bin/sh
#$ -S /bin/sh
#$ -cwd

if [ ! -d ../result/annot ];
then
    mkdir -p ../result/annot
fi

if [ ! -d ../result/22liver_nonCtr ];
then
    mkdir -p ../result/22liver_nonCtr
    if [ ! -d ../result/22liver_nonCtr/filt ];
    then
        mkdir -p ../result/22liver_nonCtr/filt
    fi
    if [ ! -d ../result/22liver_nonCtr/annot ];
    then
        mkdir -p ../result/22liver_nonCtr/annot
    fi
    if [ ! -d ../result/22liver_nonCtr/mut2sp ];
    then
        mkdir -p ../result/22liver_nonCtr/mut2sp
    fi
fi

if [ ! -d ../result/22liver_withCtr ];
then
    mkdir -p ../result/22liver_withCtr
    if [ ! -d ../result/22liver_withCtr/filt ];
    then
        mkdir -p ../result/22liver_withCtr/filt
    fi
    if [ ! -d ../result/22liver_withCtr/annot ];
    then
        mkdir -p ../result/22liver_withCtr/annot
    fi
    if [ ! -d ../result/22liver_withCtr/mut2sp ];
    then
       mkdir -p ../result/22liver_withCtr/mut2sp
    fi
fi

if [ ! -d ../result/22liver_GTEx ];
then
    mkdir -p ../result/22liver_GTEx
    if [ ! -d ../result/22liver_GTEx/filt ];
    then
        mkdir -p ../result/22liver_GTEx/filt
    fi
    if [ ! -d ../result/22liver_GTEx/annot ];
    then
        mkdir -p ../result/22liver_GTEx/annot
    fi
    if [ ! -d ../result/22liver_GTEx/mut2sp ];
    then
       mkdir -p ../result/22liver_GTEx/mut2sp
    fi
fi

:<<_COMMENT_OUT_

# 22 liver without control data
while read sample;
do
    echo "python simpleFilt.py ../data/RK/mapsplice/${sample}C.junctions.txt 5 > ../result/22liver_nonCtr/filt/${sample}C.junctions.filt.txt"
    python simpleFilt.py ../data/RK/mapsplice/${sample}C.junctions.txt 5 > ../result/22liver_nonCtr/filt/${sample}C.junctions.filt.txt

    echo " python annotSplicing.py ../result/22liver_nonCtr/filt/${sample}C.junctions.filt.txt db/refGene.bed.gz db/refExon.bed.gz > ../result/22liver_nonCtr/annot/${sample}.junction.annot.txt"
    python annotSplicing.py ../result/22liver_nonCtr/filt/${sample}C.junctions.filt.txt db/refGene.bed.gz db/refExon.bed.gz > ../result/22liver_nonCtr/annot/${sample}.junction.annot.txt
done < ../data/RK/sample_list.txt

while read sample1;
do
    while read sample2;
    do
        echo "python getSplicingMut.py ../result/22liver_nonCtr/annot/${sample1}.junction.annot.txt ../data/vcf/${sample2}.mutation.vcf.gz db/refExon.bed.gz > ../result/22liver_nonCtr/mut2sp/${sample1}-${sample2}.mutation2splicing.txt"
        python getSplicingMut.py ../result/22liver_nonCtr/annot/${sample1}.junction.annot.txt ../data/vcf/${sample2}.mutation.vcf.gz db/refExon.bed.gz > ../result/22liver_nonCtr/mut2sp/${sample1}-${sample2}.mutation2splicing.txt
    done < ../data/RK/sample_list.txt
done < ../data/RK/sample_list.txt


# 22 liver with control data
while read sample;
do
    echo "python simpleFilt.py ../data/RK/mapsplice/${sample}L.junctions.txt 2 > ../result/22liver_withCtr/filt/${sample}L.junctions.filt.txt"
    python simpleFilt.py ../data/RK/mapsplice/${sample}L.junctions.txt 2 > ../result/22liver_withCtr/filt/${sample}L.junctions.filt.txt
 done < ../data/RK/sample_list.txt

echo "cat ../result/22liver_withCtr/filt/*L.junctions.filt.txt | cut -f 1-3 | sort -k1,1 -k2,2n -k3,3n | uniq > ../result/22liver_withCtr/filt/control.merged.bed"
cat ../result/22liver_withCtr/filt/*L.junctions.filt.txt | cut -f 1-3 | sort -k1,1 -k2,2n -k3,3n | uniq > ../result/22liver_withCtr/filt/control.merged.bed

echo "bgzip -f ../result/22liver_withCtr/filt/control.merged.bed > ../result/22liver_withCtr/filt/control.merged.bed.gz"
bgzip -f ../result/22liver_withCtr/filt/control.merged.bed > ../result/22liver_withCtr/filt/control.merged.bed.gz

echo "tabix -p bed ../result/22liver_withCtr/filt/control.merged.bed.gz"
tabix -p bed ../result/22liver_withCtr/filt/control.merged.bed.gz


while read sample;
do
    echo "python filterControl.py ../data/RK/mapsplice/${sample}C.junctions.txt ../result/22liver_withCtr/filt/control.merged.bed.gz 5 > ../result/22liver_withCtr/filt/${sample}C.junctions.filt.txt"
    python filterControl.py ../data/RK/mapsplice/${sample}C.junctions.txt ../result/22liver_withCtr/filt/control.merged.bed.gz 5 > ../result/22liver_withCtr/filt/${sample}C.junctions.filt.txt

    echo "python annotSplicing.py ../result/22liver_withCtr/filt/${sample}C.junctions.filt.txt db/refGene.bed.gz db/refExon.bed.gz > ../result/22liver_withCtr/annot/${sample}.junction.annot.txt"
    python annotSplicing.py ../result/22liver_withCtr/filt/${sample}C.junctions.filt.txt db/refGene.bed.gz db/refExon.bed.gz > ../result/22liver_withCtr/annot/${sample}.junction.annot.txt

done < ../data/RK/sample_list.txt


while read sample1;
do
    while read sample2;
    do
        echo "python getSplicingMut.py ../result/22liver_withCtr/annot/${sample1}.junction.annot.txt ../data/vcf/${sample2}.mutation.vcf.gz db/refExon.bed.gz > ../result/22liver_withCtr/mut2sp/${sample1}-${sample2}.mutation2splicing.txt"
        python getSplicingMut.py ../result/22liver_withCtr/annot/${sample1}.junction.annot.txt ../data/vcf/${sample2}.mutation.vcf.gz db/refExon.bed.gz > ../result/22liver_withCtr/mut2sp/${sample1}-${sample2}.mutation2splicing.txt
    done < ../data/RK/sample_list.txt
done < ../data/RK/sample_list.txt



# 22 liver without control data
while read sample;
do
    echo "python filterGTEx.py ../data/RK/mapsplice/${sample}C.junctions.txt ../data/GTEx/Liver.junctions.bed.gz 5 2 2 > ../result/22liver_GTEx/filt/${sample}C.junctions.filt.txt"
    python filterGTEx.py ../data/RK/mapsplice/${sample}C.junctions.txt ../data/GTEx/Liver.junctions.bed.gz 5 2 2 > ../result/22liver_GTEx/filt/${sample}C.junctions.filt.txt

    echo " python annotSplicing.py ../result/22liver_GTEx/filt/${sample}C.junctions.filt.txt db/refGene.bed.gz db/refExon.bed.gz > ../result/22liver_GTEx/annot/${sample}.junction.annot.txt"
    python annotSplicing.py ../result/22liver_GTEx/filt/${sample}C.junctions.filt.txt db/refGene.bed.gz db/refExon.bed.gz > ../result/22liver_GTEx/annot/${sample}.junction.annot.txt
done < ../data/RK/sample_list.txt

_COMMENT_OUT_

while read sample1;
do
    while read sample2;
    do
        echo "python getSplicingMut.py ../result/22liver_GTEx/annot/${sample1}.junction.annot.txt ../data/RK/vcf/${sample2}.mutation.vcf.gz db/refExon.bed.gz > ../result/22liver_GTEx/mut2sp/${sample1}-${sample2}.mutation2splicing.txt"
        python getSplicingMut.py ../result/22liver_GTEx/annot/${sample1}.junction.annot.txt ../data/RK/vcf/${sample2}.mutation.vcf.gz db/refExon.bed.gz > ../result/22liver_GTEx/mut2sp/${sample1}-${sample2}.mutation2splicing.txt
    done < ../data/RK/sample_list.txt
done < ../data/RK/sample_list.txt


