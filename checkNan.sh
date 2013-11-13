#!/bin/bash

FILES=/scratch/xc406/hg19fimo/hg19gff1e3loci74/htseq_output/*_DGF_Th17_htseq_output_poissoncdf

        for file in $FILES
        do
		dirname="${file%/*}"
		dirnamep="${dirname%/*}"
		filename="${file}"
                name="${filename##*/}"
		tfname=$(echo $name |cut -d'_' -f1)
		cat $file | grep 'nan' > $dirname/temp
                if [ ! -s $dirname/temp ]; then
			echo $name
                        rm $dirnamep/poissonpbs/${tfname}_DGF_Th17_htseq_output_emppoisson.pbs
			rsync --bwlimit 30000 $file /home/xc406/hg19priors/v073113/${name}_u
			sort -k2 -n -r /home/xc406/hg19priors/v073113/${name}_u > /home/xc406/hg19priors/v073113/$name
			rm /home/xc406/hg19priors/v073113/${name}_u
                fi
        done
