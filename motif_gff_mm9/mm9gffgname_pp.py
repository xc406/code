import sys, os, csv, pp
##from multiprocessing import Pool, Process

def chunks(tflist,n):
    for i in xrange(0, len(tflist), n):
        yield tflist[i:i+n]##break tflist into chunks

def gffmod(infile, writer, motif, locidict):    
    ifile = open(infile, "rt")
    reader = csv.reader(ifile, delimiter = '\t')
    ifile.seek(0)
    for row in reader:
        if row[0] == motif:
            #print row[4]
            #row.append(row[-1])
            #nline = row[1].split('_')
            pstop = row[3]#store stop pos
            row[3] = str(int(locidict[row[1]][1])+int(row[2])-1)#start
            pval = row[-3]#store pval
            row[6] = row[4]#strand
            row[5] = pval#pval
            row[4] = str(int(locidict[row[1]][1])+int(pstop)-1)#stop        
            #nm = nline[0] + '_' + nline[1]#NM#
            motifID = row[0]#motifID
            row[0] = locidict[row[1]][0]#chr
            row[7] = 'gene_id ' + row[1] + '; sequence ' + row[-1]#store group
            row[8] = row[7]#group
            row[1] = motifID
            row[7] = '.'#frame
            row[2] = 'motif'#feature
            writer.writerows([row])

def gffgen(tfchunk, tfdict, locidict, infile):
    for n in tfchunk:
        if not tfdict[n] == '.':
       	    #if n == 'Stat3': #or n == 'IRF4' or n == 'RORC': 
            ofile = open('/scratch/xc406/mm9fimo/mm9gff1e4ud10/' + n + '.gff','w')
            writer = csv.writer(ofile, delimiter = '\t')
            if ',' in tfdict[n]:
                motifs = tfdict[n].split(',')
                #print motifs
                for m in motifs:    
                    gffmod(infile,writer,m,locidict)
            else:
                gffmod(infile,writer,tfdict[n],locidict)

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

    print "Starting pp with", job_server.get_ncpus(), "workers"
    
    infile = sys.argv[1]

    (path,fname) = os.path.split(infile)

    #ifile = open(infile,'rt')
    
    ##writer = csv.writer(ofile, delimiter = '\t')

    ##(path1,fname1) = os.path.split(infile1)
    ##(path2,fname2) = os.path.split(infile2)

    ifiletf = open('/scratch/xc406/mm9fimo/TF_Information_all_motifs.txt','rt')
    readertf = csv.reader(ifiletf, delimiter = '\t')

    ifileloci = open('/scratch/xc406/mm9fimo/mm9genesplus10kb','r')
    readerloci = csv.reader(ifileloci, delimiter = '\t')
    
    mylist1 = list([(row[3],row[6]) for row in readertf])##list of (tf,motifID)
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

    tfchunks = list(chunks(tflist, 63))
    
    #print len(tfchunks)
    #threads= []
    #global threadLock
    #threadLock = threading.Lock()
    #(self, threadID, name, tfchunks, tfdict, locidict, ifile)
    #thread1 = myThread(tfchunks, , "Thread-1", 1)
    #thread2 = myThread(2, "Thread-2", 2)
    
    ## submitting 24 jobs
    jobs = [(tfchunk,job_server.submit(gffgen,(tfchunk, tfdict, locidict, infile,),(gffmod,),("csv","os","sys",))) for tfchunk in tfchunks]
    for tfchunk, job in jobs:
	print "processing", tfchunks.index(tfchunk)
	result = job()
	if result:
	    break
    
    job_server.wait()
    job_server.print_stats

    #thread.start_new_thread( gffgen, (tfchunks[1], tfdict, locidict, ifile, ) )
    #thread.start_new_thread( gffgen, (tfchunks[2], tfdict, locidict, ifile, ) )
    #thread.start_new_thread( gffgen, (tfchunks[3], tfdict, locidict, ifile, ) )
    #for i in range(len(tfchunks)):
	#print "starting myThread" + str(i)
        #thread = myThread(i, "myThread"+str(i), tfchunks[i], tfdict, locidict, ifile)
        #thread.start()
	#threads.append(thread)
    
    #for t in threads:
	#t.join()

if __name__=='__main__':
    
##    threadLock = threading.Lock()
##    threads = []
##    for i in range(12):
##    	threads[i] = myThread(i, "Thread"+i, i)
##    	threads[i].start()
##    	threads[i].join()
    #for i in range(12):
	#p = Process(target=gffgen, args=((tfchunks[i],tfdict,locidict,ifile),))
	#p.start()
    #pool = Pool(processes = 12)
    #result = pool.apply_async(main, [])
    sys.exit(main(sys.argv))




