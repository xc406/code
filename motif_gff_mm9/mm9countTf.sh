#!/bin/bash

FILES=/home/xc406/data/mm9motifs80/pwms_all_motifs/motifoutput/*
TFFILES=/home/xc406/data/mm9motifs80/mappable/motifoutput/*.txt
for tf in $TFFILES
do 
    tfpath=${tf%/*}
    tfdir=${tf%.*}
    #echo $tfdir
    tfname="${tfdir##*/}"
    mid=$(echo $tfname |cut -d'_' -f1)
    v='_0.80'
    mname=$mid$v	
    #echo $mname
    #mkdir $tfdir #/home/xc406/data/Homo_sapiens_2013_02_20_3-54_pm/motifsPerTf/${tfname}mo
    cat /home/xc406/data/mm9motifs80/TF_Information80mm9.txt | grep ${mname} 

done

