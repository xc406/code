#!/bin/bash

FILES=/home/xc406/data/Homo_sapiens_2013_02_20_3-54_pm/motifoutput/*

for f in $FILES
do  
    filename=${f%.*}
    name="${filename##*/}"
  
    for i in $(< /home/xc406/data/Homo_sapiens_2013_02_20_3-54_pm/motifoutput/motifTh/motifTh)
    do
	if [[ "$name" == *"$i"* ]]; then
	    echo "Im here"
	    cp $f /home/xc406/data/Homo_sapiens_2013_02_20_3-54_pm/motifoutput/motifTh/
	fi   
    done
done

