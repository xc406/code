## generate upstream 10kb files using 2BitToFa and Bedtools
1. Download hg19.2bit & mm9.2bit from /goldenPath/hg19/bigZips/ & /goldenPath/mm9/bigZips/

2. Download blatScr.zip from users.soe.ucsc.edu/~kent/src/
	modify .bashrc
		'MACHTYPE=x86_64 
		export MACHTYPE 
		export PATH="$PATH":~/bin/$MACHTYPE'
	mkdir $MACHTYPE in lib	
	make
	
	Note: BLAT pre-compiled binaries for linux as 32-bit, so to compile BLAT on x86_64 Ubuntu follow:
		apt-get install build-essential
		remove the -Werror compiler flag that treats warnings as errors
			edit the /inc/common.mk
			HG_WARN_ERR = -DJK_WARN -Wall -Werrror
			to
			HG_WARN_ERR = -DJK_WARN -Wall
		makedir -p ~/bin/x86_64 
		export MACHTYPE=x86_64 (add this to $PATH)
		make
    
3. Convert 2bit to fasta file format 
	twoBitToFa mm9.2bit mm9.fa

4. Generate BED file using UCSC Tables to output upstream10kb.bed

5. Use bedtools to generate upstream10kb.fa
	bedtools getfasta -name -s -fi mm9.fa -bed mm9upstream10kb.bed -fo upstream10kb.fa
	'-s forces strand information'

## motif scan
1. Download pwms and TF_Info files by species from the Hughes database: http://cisbp2.ccbr.utoronto.ca/

2. Convert the pwms into uniprobe format matrices
	python motifs.py <path-to-the-motif-files>

3. Convert uniprobe format matrices into meme compatible format in Command-line
	uniprobe2meme -bg mm9upstream10kbbgfile ./motif_output/M*.txt > mm9motifs.meme

4. Use fimo in MEME Suite to search for alignment
	fasta-get-markov < mm9upstream10kb.fa > mm9upstream10kbbgfile
	fimo --text --bgfile upstream10kbbgfile ./mm9motifs.meme ./mm9upstream10kb.fa > mm9fimoout_date.txt 2> mm9fimoerr_date.txt

## formating fimo outputs to gff files
1. 


