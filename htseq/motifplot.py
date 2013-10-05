import HTSeq, numpy
from matplotlib import pyplot
import os, sys
import pickle

def main(argv):
    if len(argv) < 3:
        sys.stderr.write("Usage: %s bam_file gff_file window_size\n" % argv[0])
        return 1
    #if not os.path.isfile(argv[1]):
        #sys.stderr.write('Error: bam_file %r was not found!\n' % argv[1])
        #return 1
    if not os.path.isfile(argv[1]):
        sys.stderr.write('Error: gff_file %r was not found!\n' % argv[1])
        return 1
    #if not os.path.isfile(argv[2]):
        #sys.stderr.write('Error: window_size %r was not defined!\n' % argv[2])
        #return 1

    #infile_bam = sys.argv[1]
    infile_gff = sys.argv[1]

    #(bampath,bamfname) = os.path.split(infile_bam)
    
    (gffpath,gfffname) = os.path.split(infile_gff)
    (tfname,ext) = os.path.splitext(gfffname)

    #ifile_exp = open(infile_exp,'rt')
    #reader_exp = csv.reader(ifile_exp, delimiter = '\t')
    
    #bamfile = HTSeq.BAM_Reader( infile_bam )#
    bamfile = HTSeq.BAM_Reader("http://hgdownload-test.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHistone/wgEncodeBroadHistoneK562P300StdAlnRep1.bam")#wgEncodeOpenChromFaire/wgEncodeOpenChromFaireK562AlnRep1.bam" )#bamfname )#"wgEncodeUwDnaseTh17AlnRep1.bam" 
    gfffile = HTSeq.GFF_Reader( infile_gff )#"Homo_sapiens.GRCh37.70.gtf" 
    halfwinwidth = int(sys.argv[2])
    fragmentsize = 200


    ##check if a read overlaps a window
    motifpos = HTSeq.GenomicArrayOfSets( "auto", stranded=False )
		        
    for feature in gfffile:
        #print feature.iv.chrom
        #motif = HTSeq.GenomicInterval(row[0],row[3],row[4],'.')
	if feature.iv.start < halfwinwidth:
            motif = HTSeq.GenomicInterval(feature.iv.chrom, feature.iv.start-feature.iv.start, feature.iv.end+halfwinwidth, feature.iv.strand)
        else:
	    motif = HTSeq.GenomicInterval(feature.iv.chrom, feature.iv.start-halfwinwidth, feature.iv.end+halfwinwidth, feature.iv.strand)       
        motifpos [ motif ] += feature.iv.start_d_as_pos
        #print motif
	      
    #with open('motifpos.txt','w') as f:
        #pickle.dump(motifpos,f)

    #with open('motifpos.txt','r') as f:
	#motifpos = pickle.load(f) #eval(f.read())

    profile = numpy.zeros( 2*halfwinwidth, dtype="i" )
    for almnt in bamfile:
       if almnt.aligned:
          almnt.iv.length = fragmentsize
          ## take the union of all step sets that a TSS appears in
          s = set()
          #if not almnt.iv.chrom == "chrM":
          if (almnt.iv.start >= halfwinwidth) and (almnt.iv.chrom != "chrM"):
	      for step_iv, step_set in motifpos[ almnt.iv ].steps():
              #print almnt.iv.chrom, almnt.iv.strand, almnt.iv.start, almnt.iv.end, '<-- alignment'
	          print motifpos[ almnt.iv ], '<-- tsspos[alignment]' 
	          print step_iv, '<-- interval'
	          print step_set, '<- set'
	          s |= step_set
	          print s, '<-- newSet'

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


    numpy.savetxt(os.path.join(gffpath, 'p300' + tfname + 'motifplot.txt'), profile, delimiter = '\t')
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
    #numpy.savetxt(os.path.join(bampath, bamfname + tfname + 'motifplot.txt'), profile, delimiter = '\t')

# load example
#f = open( filename , 'r' )
#data = [i.strip().split( '\t' ) for i in f.xreadlines()]
#f.close()
#vector1 = [float( i[0] ) for i in data]
#vector2 = [float( i[1] ) for i in data]
#pyplot.plot( vector1 , vector2 )

if __name__=='__main__':
    sys.exit(main(sys.argv))

