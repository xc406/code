import HTSeq, numpy
from matplotlib import pyplot

bamfile = HTSeq.BAM_Reader( "wgEncodeUwDnaseTh17AlnRep1.bam" )
gtffile = HTSeq.GFF_Reader( "Homo_sapiens.GRCh37.70.gtf" )
halfwinwidth = 2000
fragmentsize = 200

##check if a read overlaps a window
tsspos = HTSeq.GenomicArrayOfSets( "auto", stranded=False )
for feature in gtffile:
   #print feature.chrom
   #break
   if feature.type == "exon" and feature.attr["exon_number"] == "1":
      p = feature.iv.start_d_as_pos
      #print p.chrom, p.pos, p 
      if not p.pos < 2000: 
          window = HTSeq.GenomicInterval( p.chrom, p.pos - halfwinwidth, p.pos + halfwinwidth, "." )
          tsspos[ window ] += p


profile = numpy.zeros( 2*halfwinwidth, dtype="i" )
for almnt in bamfile:
   if almnt.aligned:
      almnt.iv.length = fragmentsize
      ## take the union of all step sets that a TSS appears in
      s = set()
      if not almnt.iv.chrom == "chrM":
          for step_iv, step_set in tsspos[ almnt.iv ].steps():
             s |= step_set
	     #print almnt.iv.chrom, almnt.iv.strand, almnt.iv.start, p.pos, almnt.iv.end
	     #print p.chrom

      for p in s:
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
pyplot.plot( numpy.arange( -halfwinwidth, halfwinwidth ), profile )  
pyplot.savefig( 'tss.pdf')
#pyplot.show()  

# save to file
filename = 'out.txt'
f = open( filename , 'w' )
temp = numpy.arange( -halfwinwidth , halfwinwidth )
f.rwite( '\n'.join( [temp[i] +'\t'+ profile[i] for i in xrange( len( profile ) )] ) )
f.close()

# load example
#f = open( filename , 'r' )
#data = [i.strip().split( '\t' ) for i in f.xreadlines()]
#f.close()
#vector1 = [float( i[0] ) for i in data]
#vector2 = [float( i[1] ) for i in data]
#pyplot.plot( vector1 , vector2 )





















