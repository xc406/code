#!/bin/bash

FILES=/home/xc406/data/mm9phylop/*.wigFix

for f in $FILES
do

        filename=${f%.*}
        name="${filename##*/}"   
        wigToBed $f 0 ../data/mm9phylop/${name}.bed

done

