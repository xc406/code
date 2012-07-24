#!/bin/bash

FILES=/home/xc406/data/hg19gff/*.gff

for f in $FILES
do    
    if [ -s $f ]; then
    	python hg19overlap.py $f
    fi
done

FILES2=/home/xc406/data/hg19gff_gname/*.gff

for f in $FILES2
do

        filename=${f%.*}
        name="${filename##*/}"   
        sort -k1,1 -k4,4n $f > ../data/hg19gff_sorted/${name}.gff

done

FILES3=/home/xc406/data/hg19gff_sorted/*.gff

for f in $FILES3
do
    python hg19overlap2.py $f

done

FILES4=/home/xc406/data/hg19gff_final/*.gff

for f in $FILES4
do
    python hg19motif_window.py $f

done
