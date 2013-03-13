import sys, os, csv, pp
##from multiprocessing import Pool, Process
##import thread

def chunks(tflist,n):
    for i in xrange(0, len(tflist), n):
        yield tflist[i:i+n]##break tflist into chunks

def gffmod(infile, writer, motif):
    
    ifile = open(infile, "rt")
    reader = csv.reader(ifile, delimiter = '\t')
    ifile.seek(0)
    for row in reader:
	if not '_' in row[1]: 
            if row[0] == motif:
            #print row[4]
            #row.append(row[-1])
            #nline = row[1].split('_') 
                pstop = row[3]#store stop pos
                row[3] = row[2]#start
            	pval = row[-3]#store pval
            	row[6] = row[4]#strand
            	row[5] = pval#pval
            	row[4] = pstop#stop        
            #nm = nline[0] + '_' + nline[1]#NM#
            	motifID = row[0]#motifID
            	row[0] = row[1]#chr
            	row[7] = 'gene_id ' + row[1] + '_' + row[2] + '_' + str(int(pstop)-int(row[2])) + '; sequence ' + row[-1]#store group
            	row[8] = row[7]#group
            	row[1] = motifID
            	row[7] = '.'#frame
            	row[2] = 'motif'#feature
            	writer.writerows([row])
    

def gffgen(tfchunk, tfdict, infile):
    for n in tfchunk:
        if not tfdict[n] == '.':
       	    #if n in ["GATA3","RORC","TBX21","BATF","IRF4"]: 
    	    ofile = open('/scratch/xc406/hg19fimo/hg19thgff1e3loci/' + n + '.gff', 'w')
            writer = csv.writer(ofile, delimiter = '\t')
            if ',' in tfdict[n]:
                motifs = tfdict[n].split(',')
                #print motifs
                for m in motifs:
		    #ofile = open('/scratch/xc406/hg19fimo/hg19gff1e3loci/' + m + '_' + n + '.gff', 'w')
            	    #writer = csv.writer(ofile, delimiter = '\t')    
                    gffmod(infile,writer,m)
            else:
		#ofile = open('/scratch/xc406/hg19fimo/hg19gff1e3loci/' + tfdict[n] + '_' + n + '.gff', 'w')
                #writer = csv.writer(ofile, delimiter = '\t')
                gffmod(infile,writer,tfdict[n])

            ofile.close()
    
def main(argv):
    ppservers = ()

    if len(argv) < 3:
        sys.stderr.write("Usage: %s fimo_output_file #_of_workers\n" % argv[0])
        return 1
    if not os.path.isfile(argv[1]):
        sys.stderr.write('Error: fimo_output_file %r was not found!\n' % argv[1])
        return 1
    if len(sys.argv) > 2:
    	ncpus = int(sys.argv[2])
    	job_server = pp.Server(ncpus, ppservers=ppservers)
    elif len(sys.argv) == 2:
	job_server = pp.Server(ppservers=ppservers)

##    if not os.path.isfile(argv[2]):
##        sys.stderr.write('Error: TF_Info_file %r was not found!\n' % argv[2])
##        return 1
##    if not os.path.isfile(argv[3]):
##        sys.stderr.write('Error: nm#_conversion_file %r was not found!\n' % argv[3])
##        return 1

    print "Starting pp with", job_server.get_ncpus(), "workers"
    
    infile = sys.argv[1]

    (path,fname) = os.path.split(infile)

    #ifile = open(infile,'rt')
    
    ##writer = csv.writer(ofile, delimiter = '\t')

    ##(path1,fname1) = os.path.split(infile1)
    ##(path2,fname2) = os.path.split(infile2)

    ifiletf = open('/scratch/xc406/hg19fimo/TF_Information_all_motifs71.txt','rt')
    readertf = csv.reader(ifiletf, delimiter = '\t')

    #ifileloci = open('/scratch/xc406/hg19fimo/hg19ud5.bed','r')
    #readerloci = csv.reader(ifileloci, delimiter = '\t')
    
    mylist1 = []
    for row in readertf:
	if row[6] in ["STAT3","GATA3","RORC","TBX21","BATF","IRF4"]:
	    mylist1.append((row[3],row[6]))

    #mylist1 = list([(row[3],row[6]) for row in readertf])##list of (motifID,tf)
    #mylist2 = list([(row[0],int(row[1]),row[-1]) for row in readerloci])##list of (chr,start,gname)
   
    ##generate a dictionary of [tf,motifID]
    tfdict = {}
    for m,n in mylist1:
        if not n in tfdict:
            tfdict[n]= m
        else:
            tfdict[n] += ','
            tfdict[n] += m

    ##generate a dictionary of [gname,(chr,start)]
    #locidict = {}
    #for c,l,g in mylist2:
	#if not g in locidict:
	    #locidict[g]= (c,l)
    
    #print len(locidict)    
    tflist = tfdict.keys()##tf list
    #gnlist = locidict.keys()##gname list

    tfname= tfdict.items()
    ##print len(tflist) 
    ##print tfdict

    tfchunks = list(chunks(tflist, 1))
    
    ## Submit 24 jobs at once
    jobs = [(tfchunk,job_server.submit(gffgen,(tfchunk, tfdict, infile,),(gffmod,),("csv","os","sys",))) for tfchunk in tfchunks]
    for tfchunk, job in jobs:
	print "processing", tfchunks.index(tfchunk)
	result = job()
	if result:
	    break
    
    job_server.wait()
    job_server.print_stats

if __name__=='__main__':
    sys.exit(main(sys.argv))




