import HTSeq, numpy
from matplotlib import pyplot
import os, sys
import pickle

def main(argv):
    if len(argv) < 4:
        sys.stderr.write("Usage: %s bam_file gtf_file window_size\n" % argv[0])
        return 1
    if not os.path.isfile(argv[1]):
        sys.stderr.write('Error: bam_file %r was not found!\n' % argv[1])
        return 1
    if not os.path.isfile(argv[2]):
        sys.stderr.write('Error: gtf_file %r was not found!\n' % argv[2])
        return 1
    #if not os.path.isfile(argv[3]):
        #sys.stderr.write('Error: window_size %r was not defined!\n' % argv[3])
        #return 1

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
    ###tsspos = HTSeq.GenomicArrayOfSets( "auto", stranded=False )
    ###for feature in gtffile:
       #print feature.chrom
       #break
       ###if feature.type == "exon" and feature.attr["exon_number"] == "1":
          ###p = feature.iv.start_d_as_pos
          ###if not p.chrom.startswith('H'):
              ###p.chrom = 'chr' + p.chrom
          #print p.chrom, p.pos, p 
          ###if not p.pos < halfwinwidth:
	      ###if 'chr' in p.chrom: 
                  ###window = HTSeq.GenomicInterval( p.chrom, p.pos - halfwinwidth, p.pos + halfwinwidth, "." )
                  #print window
	          ###tsspos[ window ] += p
    #with open('tsspos.txt','w') as f:
        #pickle.dump(tsspos,f)

    with open('tsspos.txt','r') as f:
	tsspos = pickle.load(f) #eval(f.read())

    profile = numpy.zeros( 2*halfwinwidth, dtype="i" )
    for almnt in bamfile:
       if almnt.aligned:
          almnt.iv.length = fragmentsize
          ## take the union of all step sets that a TSS appears in
          s = set()
          #if not almnt.iv.chrom == "chrM":
          if not almnt.iv.start < halfwinwidth:
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

#fig = pyplot.figure()
#pyplot.plot( numpy.arange( -halfwinwidth, halfwinwidth ), profile )  
#pyplot.savefig( 'tss.pdf')
#pyplot.show()  

# save to file
#filename = 'out.txt'
#f = open( filename , 'w' )
#temp = numpy.arange( -halfwinwidth , halfwinwidth )
#f.write( '\n'.join( ['\t'.join([temp[i], profile[i]]) for i in numpy.arange( len( profile ) )] ) )
#f.close()
    numpy.savetxt(bampath + bamfname + 'plot.txt', profile, delimiter = '\t')

# load example
#f = open( filename , 'r' )
#data = [i.strip().split( '\t' ) for i in f.xreadlines()]
#f.close()
#vector1 = [float( i[0] ) for i in data]
#vector2 = [float( i[1] ) for i in data]
#pyplot.plot( vector1 , vector2 )

if __name__=='__main__':
    sys.exit(main(sys.argv))



















