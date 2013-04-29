import HTSeq, numpy
from matplotlib import pyplot
import os, sys, pp
import pickle

def chunks(bamfile,n):
    for i in xrange(0, len(bamfile), n):
        yield bamfile[i:i+n]##break bamlist into chunks

def tssposdict(gtffile,halfwinwidth):
    ##check if a read overlaps a window
    tsspos = HTSeq.GenomicArrayOfSets( "auto", stranded=False )
    for feature in gtffile:
       #print feature.chrom
       #break
       if feature.type == "exon" and feature.attr["exon_number"] == "1":
          p = feature.iv.start_d_as_pos
          if not p.chrom.startswith('H'):
              p.chrom = 'chr' + p.chrom
          #print p.chrom, p.pos, p 
          if not p.pos < halfwinwidth:
              if 'chr' in p.chrom: 
                  window = HTSeq.GenomicInterval( p.chrom, p.pos - halfwinwidth, p.pos + halfwinwidth, "." )
                  #print window
                  tsspos[ window ] += p
    return tsspos

def profilegen(bamfile,halfwinwidth,fragmentsize):
    profile = numpy.zeros( 2*halfwinwidth, dtype="i" )
    for almnt in bamfile:
       if almnt.aligned:
          almnt.iv.length = fragmentsize
          ## take the union of all step sets that a TSS appears in
          s = set()
          #if not almnt.iv.chrom == "chrM":
          if (almnt.iv.start >= halfwinwidth) and (almnt.iv.chrom != "chrM"):
            for step_iv, step_set in tsspos[ almnt.iv ].steps():
                 print almnt.iv.chrom, almnt.iv.strand, almnt.iv.start, almnt.iv.end
                 #if almnt.iv.strand == '+':
                 #if not almnt.iv.start < halfwinwidth: 
                 s |= step_set
                 #if almnt.iv.strand == '-':
                     #if not almnt.iv.start > p.pos + halfwinwidth:
                        #s |= step_set 

                 #print almnt.iv.chrom, almnt.iv.strand, almnt.iv.start, p.pos, almnt.iv.end
                 #print p.chrom

          for p in s:
             #print p
             #if not almnt.iv.chrom == "chrM":
             if p.strand == "+":
                start_in_window = almnt.iv.start - p.pos + halfwinwidth
                end_in_window   = almnt.iv.end   - p.pos + halfwinwidth
             else:
                start_in_window = p.pos + halfwinwidth - almnt.iv.end
                end_in_window   = p.pos + halfwinwidth - almnt.iv.start
             start_in_window = max( start_in_window, 0 )
             end_in_window = min( end_in_window, 2*halfwinwidth )
             profile[ start_in_window : end_in_window ] += 1
    return profile	

def main(argv):
    ppservers = ()
    if len(argv) < 4:
        sys.stderr.write("Usage: %s bam_file gtf_file window_size #_of_workers\n" % argv[0])
        return 1
    if not os.path.isfile(argv[1]):
        sys.stderr.write('Error: bam_file %r was not found!\n' % argv[1])
        return 1
    if not os.path.isfile(argv[2]):
        sys.stderr.write('Error: gtf_file %r was not found!\n' % argv[2])
        return 1
    if len(sys.argv) == 5:
        ncpus = int(sys.argv[4])
        job_server = pp.Server(ncpus, ppservers=ppservers)
    elif len(sys.argv) == 4:
        job_server = pp.Server(ppservers=ppservers)
    #if not os.path.isfile(argv[3]):
        #sys.stderr.write('Error: window_size %r was not defined!\n' % argv[3])
        #return 1

    print "Starting pp with", job_server.get_ncpus(), "workers"

    infile_bam = sys.argv[1]
    infile_gtf = sys.argv[2]

    (bampath,bamfname) = os.path.split(infile_bam)
    (gtfpath,gtffname) = os.path.split(infile_gtf)
    
    #ifile_exp = open(infile_exp,'rt')
    #reader_exp = csv.reader(ifile_exp, delimiter = '\t')
    
    bamfile = HTSeq.BAM_Reader( bamfname )#"wgEncodeUwDnaseTh17AlnRep1.bam" 
    gtffile = HTSeq.GFF_Reader( gtffname )#"Homo_sapiens.GRCh37.70.gtf" 
    halfwinwidth = int(sys.argv[3])
    fragmentsize = 200

    ##check if a read overlaps a window
    #tsspos = tssposdict(gtffile,halfwinwidth)
    #with open('tsspos.txt','w') as f:
        #pickle.dump(tsspos,f)

    with open('tsspos.txt','r') as f:
	tsspos = pickle.load(f) #eval(f.read())

    #profile = profilegen(bamfile,halfwinwidth,fragmentsize)

    #fig = pyplot.figure()
    #pyplot.plot( numpy.arange( -halfwinwidth, halfwinwidth ), profile )  
    #pyplot.savefig( 'tss.pdf')
    #pyplot.show()  

    #numpy.savetxt(bampath + bamfname + 'plot.txt', profile, delimiter = '\t')

    bamchunks = list(chunks(bamfile, 1000000))

    ## Submit 24 jobs at once
    jobs = [(tfchunk,job_server.submit(profilegen,(bamchunk, halfwinwidth, fragmentsize,),(tssposdict,),("HTSeq","numpy","os","sys",))) for bamchunk in bamchunks]
    presult = numpy.zeros( 2*halfwinwidth, dtype="i" )
    for bamchunk, job in jobs:
        print "processing", bamchunks.index(bamchunk)
        result = job()
        if result:
	    result = numpy.add(result,presult)
	    presult = result
            break

    numpy.savetxt(bampath + bamfname + 'plot.txt', result, delimiter = '\t')
    job_server.wait()
    job_server.print_stats

if __name__=='__main__':
    sys.exit(main(sys.argv))



















