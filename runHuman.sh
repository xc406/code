#!/bin/bash

cd hg19/
mkdir uniprobe
cd uniprobe
python ../../scripts/convertPWMStoUniprobe.py ../pwms/
cd ..
uniprobe2meme uniprobe/M*.txt > matrix.meme
python ../scripts/processUpstream.py upstream1000.fa
mv nmupstream1000.fa processed_upstream1000.fa
fasta-get-markov <upstream1000.fa> upstream1000.background
fimo --bgfile upstream1000.background --text ./matrix.meme ./processed_upstream1000.fa > motif_locations.fimo
mv fimoutall fimo_formatted_all
python ../scripts/format_fimo.py motif_locations.fimo TF_Information.txt hg19_refseq_May_1_2011.txt
cat hg19 | awk 'BEGIN{FS="\t";} {print $2,"\t",$4,"\t",$9} END{}'> hg.priors
