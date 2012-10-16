Input Files:
hg19_refseq_May_1_2011.txt  pwms/  TF_Information.txt  upstream1000.fa


1. Download pwms and TF_Info files by species from the Hughes database: http://cisbp2.ccbr.utoronto.ca/

2. Install Meme suite and export path

3. Convert the pwms into uniprobe format matrixes
	python motifs.py <path-to-the-motif-files>

4. Convert uniprobe format matrixes into meme compatible format in Command-line
	uniprobe2meme M*.txt > M.meme

5. Process upstream1000.fa file from ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/
    [Sequences 1000 bases upstream of annotated
    transcription starts of RefSeq genes with annotated 5' UTRs.
    This file is updated weekly so it might be slightly out of sync with
    the RefSeq data which is updated daily for most assemblies.]

	python nm.py <upstream1000.fa>
  
  
Pre-6 Generate background file
	fasta-get-markov < upstream1000.fa > upstream1000bfile

6. Do a motif search with FIMO in Command-line
	fimo --bfile upstream1000bfile --text ./M.meme ./upstream1000.fa > nmfimo.txt
	
7. Formatting prior
	python format.py <fimo_output_file> <tf_info_file> <nm#_coversion_file>
	
	A. write all columns in nmfimo.txt to fimoutall (duplicate the Motif_ID and NM# columns)
	B. substitute 2nd (Motif_ID) and 4th (NM#) column of fimoutall with informations from TF_Info.txt (TF_Names) and refseq.txt (Gene_Names) respectively
	note: both Motif_ID/TF_Names and NM#/Gene_Names are many-to-many relationship. The substitution step replaces each Motif_ID or NM# with all its corresponding TF_Names or Gene_Names separated by commas
	*C. generate a list of TFs-- tflist.py (for ENCODE Chip-seq)
*8. Verify the tf-gene pairs in JASPAR-- jaspar.py
	A. choose a gene x of interest
	B. run jaspar.py with gene x (all TFs that target gene x)
	C. search the JASPAR motif database against the upstream1000.fa of gene x
	D. compare the results in B and C 

*not part of the prior

To get a list of tfs

