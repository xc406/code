import sys, os, csv, pp
##from multiprocessing import Pool, Process
##import thread

def chunks(tflist,n):
    for i in xrange(0, len(tflist), n):
        yield tflist[i:i+n]##break tflist into chunks

def gffmod(infile, writer, motif, locidict):
    
    ifile = open(infile, "rt")
    reader = csv.reader(ifile, delimiter = '\t')
    ifile.seek(0)
    for row in reader:
	#print row
    ##fimo_infile[motif_ID,gname,local_start,local_end,strand,score,pval,matched_seq]
    ##gff_outfile[chr,motif_ID,'motif',start,end,pval,strand,.,group(gene_id,gene_name,matched_seq)]	
	if not '_' in row[1]: 
	    #print row
	    if (row[0] == motif) and (row[1] in locidict):
                newrow=[locidict[row[1]][0],row[0],'motif',str(int(locidict[row[1]][1])+int(row[2])-1),str(int(locidict[row[1]][1])+int(row[3])-1),row[6],row[4],'.','gene_id '+row[1]+'_'+locidict[row[1]][0]+'_'+str(int(locidict[row[1]][1])+int(row[2])-1)+'_'+str(int(row[3])-int(row[2])+1)+'; sequence '+row[-1]]  
      	        #print newrow
		writer.writerows([newrow])
    

def gffgen(tfchunk, tfdict, locidict, infile):
    for n in tfchunk:
        if not tfdict[n] == '.':
       	    #if n in ["GATA3","RORC","TBX21","BATF","IRF4"]: 
    	    ofile = open('/scratch/xc406/hg19fimo/hg19gff1e3loci80ud10/' + n + '.gff', 'w')
            writer = csv.writer(ofile, delimiter = '\t')
            if ',' in tfdict[n]:
                motifs = tfdict[n].split(',')
                #print motifs
                for m in motifs:
		    #ofile = open('/scratch/xc406/hg19fimo/hg19gff1e3loci/' + m + '_' + n + '.gff', 'w')
            	    #writer = csv.writer(ofile, delimiter = '\t')    
                    gffmod(infile,writer,m,locidict)
            else:
		#ofile = open('/scratch/xc406/hg19fimo/hg19gff1e3loci/' + tfdict[n] + '_' + n + '.gff', 'w')
                #writer = csv.writer(ofile, delimiter = '\t')
                gffmod(infile,writer,tfdict[n],locidict)

            ofile.close()
    
def main(argv):
    ppservers = ()

    if len(argv) < 5:
        sys.stderr.write("Usage: %s fimo_output_file #_of_workers TF_info_file loci_bed_file\n" % argv[0])
        return 1
    if not os.path.isfile(argv[1]):
        sys.stderr.write('Error: fimo_output_file %r was not found!\n' % argv[1])
        return 1
    if len(sys.argv) > 2:
    	ncpus = int(sys.argv[2])
    	job_server = pp.Server(ncpus, ppservers=ppservers)
    elif len(sys.argv) == 2:
	job_server = pp.Server(ppservers=ppservers)
    if not os.path.isfile(argv[3]):
        sys.stderr.write('Error: TF_info_file %r was not found!\n' % argv[3])
        return 1
    if not os.path.isfile(argv[4]):
        sys.stderr.write('Error: loci_info_file %r was not found!\n' % argv[4])
        return 1

    print "Starting pp with", job_server.get_ncpus(), "workers"
    
    infile = sys.argv[1]
    infiletf = sys.argv[3]
    infileloci = sys.argv[4]
    (path,fname) = os.path.split(infile)

    #ifile = open(infile,'rt')
    ##writer = csv.writer(ofile, delimiter = '\t')

    ##(path1,fname1) = os.path.split(infile1)
    ##(path2,fname2) = os.path.split(infile2)

    #ifiletf = open('/scratch/xc406/hg19fimo/TF_Information80hg19.txt','rt')
    ifiletf = open(infiletf,'rt')
    readertf = csv.reader(ifiletf, delimiter = '\t')

    #ifileloci = open('/scratch/xc406/hg19fimo/hg19se10.bed','r')
    ifileloci = open(infileloci,'r')
    readerloci = csv.reader(ifileloci, delimiter = '\t')
    
    ##mylist1 = []
    ##for row in readertf:
    	#if row[6] in ["STAT3","GATA3","RORC","TBX21","BATF","IRF4"]:
    ##	mylist1.append((row[3],row[6]))

    mylist1 = list([(row[3],row[6]) for row in readertf])##list of (motifID,tf)
    mylist2 = list([(row[0],int(row[1]),row[-1]) for row in readerloci])##list of (chr,start,gname)
   
    ##generate a dictionary of [tf,motifID]
    tfdict = {}
    for m,n in mylist1:
        if not n in tfdict:
            tfdict[n]= m
        else:
            tfdict[n] += ','
            tfdict[n] += m

    ##generate a dictionary of [gname,(chr,start)]
    locidict = {}
    for c,l,g in mylist2:
	if not g in locidict:
	    locidict[g]= (c,l)
    
    #print len(locidict)    
    tflist = tfdict.keys()##tf list
    gnlist = locidict.keys()##gname list

    tfname= tfdict.items()
    ##print len(tflist) 
    ##print tfdict

    tfchunks = list(chunks(tflist, 60))
    
    ## Submit 24 jobs at once
    jobs = [(tfchunk,job_server.submit(gffgen,(tfchunk, tfdict, locidict, infile),(gffmod,),("csv","os","sys",))) for tfchunk in tfchunks]
    for tfchunk, job in jobs:
	print "processing", tfchunks.index(tfchunk)
	result = job()
	if result:
	    break
    
    job_server.wait()
    job_server.print_stats

if __name__=='__main__':
    sys.exit(main(sys.argv))




