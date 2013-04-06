#!/bin/bash

FILES=/home/xc406/data/Homo_sapiens_2013_02_20_3-54_pm/motifoutput/*
TFFILES=/home/xc406/data/Homo_sapiens_2013_02_20_3-54_pm/motifsPerTf/*.txt
for tf in $TFFILES
do 
    tfpath=${tf%/*}
    tfdir=${tf%.*}
    #echo $tfdir
    tfname="${tfdir##*/}"
    #echo $tfname
    mkdir $tfdir #/home/xc406/data/Homo_sapiens_2013_02_20_3-54_pm/motifsPerTf/${tfname}mo
    for f in $FILES
    do  
        filename=${f%.*}
        name="${filename##*/}"
  
        for i in $(< ${tf})
        do
	    if [[ "$name" == *"$i"* ]]; then
	        echo "Im here"
	        cp $f ${tfdir}/
	    fi   
        done
    done
    uniprobe2meme -bg ~/data/hg19bgfile $tfdir/*.txt > $tfpath/hg19${tfname}.meme
done

