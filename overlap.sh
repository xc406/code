#!/bin/bash

#FILES=/home/xc406/data/mm9gff_corrected/*.gff

#for f in $FILES
#do    
#    if [ -s $f ]; then
#    	python overlap.py $f
#    fi
#done

FILES2=/home/xc406/data/mm9gff_corrected_gname/*.gff

for f in $FILES2
do

        filename=${f%.*}
        name="${filename##*/}"   
        sort -k1,1 -k4,4n $f > ../data/mm9gff_corrected_sorted/${name}.gff

done

FILES3=/home/xc406/data/mm9gff_corrected_sorted/*.gff

for f in $FILES3
do
    python overlap_2.py $f

done

FILES4=/home/xc406/data/mm9gff_corrected_final/*.gff

for f in $FILES4
do
    python motif_window.py $f

done
