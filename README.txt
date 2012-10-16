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
	note: add following commands in .bashrc
		'export PERL5LIB=/home/xc406/tools/meme/lib/perl:$PERL5LIB'	

4. Use fimo in MEME Suite to search for alignment
	fasta-get-markov < mm9upstream10kb.fa > mm9upstream10kbbgfile
	fimo --text --bgfile mm9upstream10kbbgfile --output-pthresh 1e-3 ./mm9motifs.meme ./mm9upstream10kb.fa > mm9fimoout_date.txt 2> mm9fimoerr_date.txt
	note: the default cutoff pval is 1e-4 without the --output-pthresh option
	***--psp

## formating fimo outputs to gff files
1. format fimo_output_file into gff (time-consuming)
	python mm9gff.py fimo_output_file

2. substitute/add Hugo_gene_names next to NM# (time-consuming) and take out overlaps 
	python overlap.py gff_file
	sort the gff_files
	python overlap_2.py gff_file
	--this is done with a shell script overlap.sh to first check whether the tf file is empty

3. filter by DHS reads (htseq-count_output_files)
	python filter htseq-count_output_file gff_file fimo_stderr_file
	note: unnecessary

3'. alternatively, use macs to perform peak calling with an arbitrary cutoff p val (say 1e-3) on downloaded bam files
	macs14 -t wgEncodeUwDnaseMelC57bl6MAdult8wksAlnRep1.bam -f BAM -g mm -p 1e-3 -n mel_dhs1
		note: this step can sometimes be substituted with the broadpeak/narrowpeak bed files from ENCODE
	bed files will be processed by 
		bedtools intersect -a gff_file -b macs_output_file -wa -wb > intersect_file
	process the intersect files with
		python bed_macs_two.py intersect_file

4. calculate aupr with combine.all.R
	
