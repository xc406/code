#!/bin/bash

FILES=/scratch/xc406/hg19fimo/hg19gff1e3loci/*.gff

for i in $FILES
do      
        filename="${i%.*}"
	name="${filename##*/}"
	dirname="${filename%/*}"
	cd ${dirname}/pbs
    	echo "#PBS -V
#PBS -r n
#PBS -d $dirname
#PBS -m e 
#PBS -M xc406@nyu.edu 
#PBS -S /bin/bash 
#PBS -l nodes=1:ppn=1,walltime=12:00:00,mem=4gb
#PBS -N hg19${name}gff
#PBS -q default
#PBS -l qos=quarter

source /etc/profile.d/modules.sh;
module load python/gnu/2.7.3
module load bash/gnu/4.2
	
	if [ -s $name.gff ]; then
	    #sort -k1,1 -k4n,4 $name.gff > ${dirname}/hg19gff_sorted/${name}.gff

	    python /home/xc406/code/motif_gff_hg19/hg19rmvOverlap.py $dirname/${name}.gff

            #rm $dirname/hg19gff_sorted/$name.gff
	
	    python /home/xc406/code/motif_gff_hg19/hg19motif_window.py $dirname/hg19gff_final/${name}.gff 100

	    python /home/xc406/code/motif_gff_hg19/hg19motif_windowLR.py $dirname/hg19gff_final/${name}.gff

	fi
" > $name.pbs
	cd ..
done

