import os, sys, csv, HTSeq, numpy
import pickle

def main(argv):
    if len(argv) < 5:
        sys.stderr.write("Usage: %s gtf_file gff_file left_window_size right_window_size\n" % argv[0])
        return 1
    if not os.path.isfile(argv[1]):
        sys.stderr.write('Error: gtf_file %r was not found!\n' % argv[1])
        return 1
    if not os.path.isfile(argv[2]):
        sys.stderr.write('Error: gff_file %r was not found!\n' % argv[2])
        return 1
    #if not os.path.isfile(argv[3]):
        #sys.stderr.write('Error: window_size %r was not defined!\n' % argv[3])
        #return 1

    infile_gtf = sys.argv[1]
    infile_gff = sys.argv[2]

    (gffpath,gfffname) = os.path.split(infile_gff)
    (tfname,ext) = os.path.splitext(gfffname) 
    #(gtfpath,gtffname) = os.path.split(infile_gtf)

    with open(infile_gff,'rt') as ifile_gff:
        reader_gff = csv.reader(ifile_gff, delimiter = '\t')

    ofile = open(os.path.join(gffpath,tfname + '_tss' + ext), 'w')
    writer = csv.writer(ofile, delimiter = '\t')

    #bamfile = HTSeq.BAM_Reader( "http://hgdownload-test.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeOpenChromFaire/wgEncodeOpenChromFaireK562AlnRep1.bam" )#bamfname )#"wgEncodeUwDnaseTh17AlnRep1.bam"
    gfffile = HTSeq.GFF_Reader( infile_gff )#"Homo_sapiens.GRCh37.70.gtf"
    gtffile = HTSeq.GFF_Reader( infile_gtf )#"Homo_sapiens.GRCh37.70.gtf" )
    halfwinwidth_l = int(sys.argv[3])
    halfwinwidth_r = int(sys.argv[4])
    #fragmentsize = 200

    ##check if a read overlaps a window
    tsspos = HTSeq.GenomicArrayOfSets( "auto", stranded=False )
    for feature in gtffile:
       if feature.type == "exon" and feature.attr["exon_number"] == "1":
          p = feature.iv.start_d_as_pos
	  #print p.chrom
          if not p.chrom.startswith('G'):
              p.chrom = 'chr' + p.chrom
          #print p.chrom, p.pos, p
	  if 'chr' in p.chrom:
              if p.pos < halfwinwidth_l:
	          window = HTSeq.GenomicInterval( p.chrom, p.pos - p.pos, p.pos + halfwinwidth_r, ".") 
              else:
                  window = HTSeq.GenomicInterval( p.chrom, p.pos - halfwinwidth_l, p.pos + halfwinwidth_r, "." )
              tsspos[ window ] += p
	      #print p, '<--- GenmoicPosition'
	      #print window, '<--- GenomicInterval'
	      #print tsspos[window], '<--- GenomicArrayOfSets'
    with open('mm9tsspos' + str(halfwinwidth_l) + '_' + str(halfwinwidth_r) + '.txt','w') as f:
        pickle.dump(tsspos,f)

    #with open('mm9tsspos' + str(halfwinwidth_l) + '_' + str(halfwinwidth_r) + '.txt','r') as f:
        #tsspos = pickle.load(f) #eval(f.read())

    #for row in reader_gff:
    for feature in gfffile:
	#print feature.iv.chrom
        #motif = HTSeq.GenomicInterval(row[0],row[3],row[4],'.')
	motif = HTSeq.GenomicInterval(feature.iv.chrom, feature.iv.start, feature.iv.end, feature.iv.strand)
        #print motif
        s = set()
        #if not almnt.iv.chrom == "chrM":
        #if (feature.iv.start >= halfwinwidth) and (feature.iv.chrom != "chrM"):
        for step_iv, step_set in tsspos[ motif ].steps():
                 #print almnt.iv.chrom, almnt.iv.strand, almnt.iv.start, almnt.iv.end
            #print list(tsspos[ motif ].steps())
	    if step_set == set([]):
		    #print list(feature.attr.items())
		feature.attr["sequence"] = feature.attr["sequence"].split('\r')[0]
		writer.writerows([[feature.iv.chrom, feature.source, feature.type, feature.iv.start, feature.iv.end, feature.score, feature.iv.strand, '.', ';'.join([ k + ' ' + v for k,v in feature.attr.iteritems() ])  ]])
                 #s |= step_set

    ofile.close()

if __name__=='__main__':
    sys.exit(main(sys.argv))
