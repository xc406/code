#!/bin/bash

FILES=/scratch/xc406/1000genomes/*.tbi

for f in $FILES
do    
    if [ -s $f ]; then
        echo $f
        filename="${f%.*}"
	echo $filename
	name="${filename##*/}"
	if [ -s $filename ]
	then
		echo "$name downloaded"
	else
		echo "downloading $name"
    		wget -P /scratch/xc406/1000genomes/ ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20110521/${name}
	fi
    fi
done

#FILES2=/scratch/xc406/hg19fimo/hg19gff1e3ud5/hg19gff_sorted/*.gff

#for f in $FILES2
#do

#    python /home/xc406/code/motif_gff_hg19/hg19overlap2.py $f
 
#done

#FILES3=/scratch/xc406/hg19fimo/hg19gff1e3ud5/hg19gff_final/*.gff

#for f in $FILES3
#do

#    python /home/xc406/code/motif_gff_hg19/hg19motif_windowL.py $f
#    python /home/xc406/code/motif_gff_hg19/hg19motif_windowR.py $f

#done
