##Finding enriched tf in th17 tf clusters

1. Run chipseq macs files on scorch.r to get mtls file and cluster assignment file
	python getcluster.py MTLS_File Cluster_Assignments

2. Take the mtls in cluster1 and write to fasta files with BSgenome.R
   expand each mtl by 100bp upstream and downstream respectively

3. Generate a random background 10 times the size of cluster1 
	python cluster_background.py cluster_file

4. Run motif search tool Fimo on cluster1.fa and background.fa respectively
	fimo --text --bgfile mm9bgfile --output-pthresh 1e-3 ./mm9motifs.meme ./cluster1.fa > mm9cluster1fimoout_date.txt 2> mm9cluster1fimoerr_date.txt
	fimo --text --bgfile mm9bgfile --output-pthresh 1e-3 ./mm9motifs.meme ./background.fa > mm9backgroundfimoout_date.txt 2> mm9backgroundfimoerr_date.txt

5. Run hypergeometric test on cluster vs background mtl counts and generate a file containing motif ids and p-values from phyper with motifphyper.R

6. Convert motif ids to tf names
	python motifidtogene.py motif_list TF_Info_file

7. Take the overlap between list of tfs in the Hughes database and list of tfs in the th17 core network and run Kolmogorov-Smirnov test

8. 	 
   
