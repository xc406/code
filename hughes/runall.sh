#!/usr/bin/bash

python motifs.py <path-to-the-motif-files>
python nm.py <upstream1000.fa>

##Install Meme suite
uniprobe2meme M*.txt > M.meme
fasta-get-markov < upstream1000.fa > upstream1000bfile
fimo --bfile upstream1000bfile --text ./M.meme ./upstream1000.fa > nmfimo.txt

python format.py <fimo_output_file> <tf_info_file> <nm#_coversion_file>

python tflist.py
python jaspar.py
