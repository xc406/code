import sys, os, csv, pp, io, bz2, mmap
##from multiprocessing import Pool, Process
##import thread

def chunks(tflist,n):
    for i in xrange(0, len(tflist), n):
        yield tflist[i:i+n]##break tflist into chunks

def gffmod(infile, ofile, motif):
    #ifile = open(infile, "rt")
    #reader = mmap.mmap(ifile.fileno(),0,prot=mmap.PROT_READ)
    #reader = csv.reader(ifile, delimiter = '\t')
    with bz2.BZ2File(infile,'r') as reader:
        #reader.seek(0)
        flag = True
        while flag:
	    row = reader.readline().split('\n')[0].split('\t')
	    if not row is None:
		if len(row) > 1:
	        #print row
                ##fimo_infile[motif_ID,chr,start,end,strand,score,pval,matched_seq]
                ##gff_outfile[chr,motif_ID,'motif',start,end,pval,strand,.,group(gene_id,gene_name,matched_seq)]	
	            if 'chr' in row[1]:
	                if row[0] == motif:
                            newrow=[row[1],row[0],'motif',row[2],row[3],row[6],row[4],'.','gene_id '+row[1]+'_'+row[2]+'_'+row[3]+'_'+str(int(row[3])-int(row[2])-1) +'_'+row[4]+'; sequence '+row[-1]]  
		            #print newrow
	                    ofile.writelines('\t'.join(newrow)+'\n')
      	        else:
		    flag = False
	    else:
	        flag = False

def gffgen(tfchunk, tfdict, infile, shortname):
    for n in tfchunk:
        if not tfdict[n] == '.':
       	    #if n in ["GATA3","RORC","TBX21","BATF","IRF4"]: 
    	    with bz2.BZ2File('/scratch/xc406/hg19fimo/hg19gff1e3wg80/' + n + '_' + shortname + '.gff.bz2', 'w') as ofile:
                #writer = csv.writer(ofile, delimiter = '\t')
                if ',' in tfdict[n]:
                    motifs = tfdict[n].split(',')
                    #print motifs
                    for m in motifs:
                        gffmod(infile,ofile,m)
                else:
                    gffmod(infile,ofile,tfdict[n])

            #ofile.close()
    
def main(argv):
    ppservers = ()

    if len(argv) < 4:
        sys.stderr.write("Usage: %s fimo_output_file #_of_workers TF_info_file\n" % argv[0])
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
    #if not os.path.isfile(argv[4]):
        #sys.stderr.write('Error: loci_bed_file %r was not found!\n' % argv[4])
        #return 1

    print "Starting pp with", job_server.get_ncpus(), "workers"
    
    infile = sys.argv[1]
    infiletf = sys.argv[3]
    #infileloci = sys.argv[4]

    #reader = bz2.BZ2File(infile,'r')
    (path,fname) = os.path.split(infile)
    (shortname,ext) = os.path.splitext(fname)
    #ifile = open(infile,'rt')
    ##writer = csv.writer(ofile, delimiter = '\t')

    ##(path1,fname1) = os.path.split(infile1)
    ##(path2,fname2) = os.path.split(infile2)

    #ifiletf = open('/scratch/xc406/hg19fimo/TF_Information80hg19.txt','rt')
    #ifiletf = open(infiletf,'rt')
    #readertf = csv.reader(ifiletf, delimiter = '\t')

    #ifileloci = open('/scratch/xc406/hg19fimo/hg19se10.bed','r')
    #ifileloci = open(infileloci,'r')
    #readerloci = csv.reader(ifileloci, delimiter = '\t')
    
    ##mylist1 = []
    ##for row in readertf:
    	#if row[6] in ["STAT3","GATA3","RORC","TBX21","BATF","IRF4"]:
    ##	mylist1.append((row[3],row[6]))
    with open(infiletf,'rt') as ifiletf:
        readertf = csv.reader(ifiletf, delimiter = '\t')
        mylist1 = list([(row[3],row[6]) for row in readertf])##list of (motifID,tf)
    #mylist2 = list([(row[0],int(row[1]),int(row[2]),row[-1]) for row in readerloci])##list of (chr,start,gname)
   
    ##generate a dictionary of [tf,motifID]
    tfdict = {}
    for m,n in mylist1:
        if not n in tfdict:
            tfdict[n]= m
        else:
            tfdict[n] += ','
            tfdict[n] += m

    ##generate a dictionary of [gname,(chr,start,end)]
    #locidict = {}
    #for c,s,e,g in mylist2:
	#if not g in locidict:
	    #locidict[g]= (c,s,e)
    
    #print len(locidict)    
    tflist = tfdict.keys()##tf list
    #gnlist = locidict.keys()##gname list

    tfname= tfdict.items()
    ##print len(tflist) 
    ##print tfdict

    tfchunks = list(chunks(tflist, 60))
    
    ## Submit 24 jobs at once
    jobs = [(tfchunk,job_server.submit(gffgen,(tfchunk, tfdict, infile, shortname),(gffmod,),("csv","os","sys","io","bz2"))) for tfchunk in tfchunks]
    for tfchunk, job in jobs:
	print "processing", tfchunks.index(tfchunk)
	result = job()
	if result:
	    break
    
    job_server.wait()
    job_server.print_stats

if __name__=='__main__':
    sys.exit(main(sys.argv))




