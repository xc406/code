import HTSeq, numpy
from matplotlib import pyplot
import os, sys, csv#, pp
import pickle
from collections import defaultdict

def defFeature(firewindow,feature,halfwinwidth):
    chrom = feature.split('_')[0]
    start = int(feature.split('_')[1])
    end = int(feature.split('_')[2])
    summit = abs(start-end)/2 + min(start, end)
    if summit < halfwinwidth:
        #fire = HTSeq.GenomicInterval(feature.iv.chrom, 0, summit+halfwinwidth, feature.iv.strand)
        firewindownew = HTSeq.GenomicArray("auto", stranded = True)
        for i in xrange(0,len(range(0, summit+halfwinwidth, (summit + halfwinwidth)/100))):
            firewindowiv = HTSeq.GenomicInterval(chrom,range(0, summit+halfwinwidth, (summit + halfwinwidth)/100)[i],range(0, summit+halfwinwidth, (summit + halfwinwidth)/100)[i]+(summit + halfwinwidth)/100,'-')
                #firewindowset = HTSeq.GenomicArrayOfSets("auto",stranded=False)
	    for iv, val in firewindow[firewindowiv].steps():
		firewindownew [ firewindowiv ] += val
            #firewindow [ firewindowiv ] += firewindow[ firewindowiv ]
	#for iv, val in firewindownew[firewindowiv].steps():
	    #print iv,val
	return firewindownew#firewindow[firewindowiv].steps()
    else:
        #fire = HTSeq.GenomicInterval(feature.iv.chrom, summit-halfwinwidth, summit+halfwinwidth, feature.iv.strand)
        firewindownew = HTSeq.GenomicArray("auto", stranded = True)
        for j in xrange(0,len(range(summit-halfwinwidth, summit+halfwinwidth, halfwinwidth/50))):
            firewindowiv = HTSeq.GenomicInterval(chrom,range(summit-halfwinwidth, summit+halfwinwidth, halfwinwidth/50)[j],range(summit-halfwinwidth, summit+halfwinwidth, halfwinwidth/50)[j]+halfwinwidth/50,'-')
            #print firewindowiv
	    for iv, val in firewindow[firewindowiv].steps():
		firewindownew [firewindowiv] += val
                #fire.append(firewindow)
            #print list(firewindownew[firewindowiv].steps())
	#for iv, val in firewindownew[firewindowiv].steps():
	    #print iv,val
	#print firewindow[firewindowiv]
        return firewindownew#firewindow[firewindowiv].steps()

def count(firechunk, bamfile,halfwinwidth):
    for fire in firechunk:
        #print type(firechunk[ fire ])
	#print fire
        #for firewindow in xrange(0,len(firechunk[ fire ])):
        for almnt in bamfile [fire]:
            if almnt.aligned:
                if almnt.iv.start >= halfwinwidth:
		#print almnt.iv.chrom, almnt.iv.strand, almnt.iv.start, almnt.iv.end
		    for firewindow in firechunk[fire]:#xrange(0,len(firechunk[fire])):
	                firewindow[almnt.iv] += 1
    return firechunk

def main(argv):
    if len(argv) < 3:
        sys.stderr.write("Usage: %s bam_file gff_file window_size\n" % argv[0])
        return 1
    if not os.path.isfile(argv[1]):
        sys.stderr.write('Error: bam_file %r was not found!\n' % argv[1])
        return 1
    if not os.path.isfile(argv[2]):
        sys.stderr.write('Error: gff_file %r was not found!\n' % argv[2])
        return 1

    infile_bam = sys.argv[1]
    infile_gff = sys.argv[2]

    (bampath,bamfname) = os.path.split(infile_bam)
    
    (gffpath,gfffname) = os.path.split(infile_gff)
    (tfname,ext) = os.path.splitext(gfffname)

    ofile = open(os.path.join(gffpath,'groseqplotRep3n'),'wt')
    writer = csv.writer(ofile, delimiter = '\t')

    bamfile = HTSeq.BAM_Reader( infile_bam )#
    #bamfile = HTSeq.BAM_Reader("wgEncodeUwDnaseTh17AlnRep1.bam") 
    gfffile = HTSeq.GFF_Reader( infile_gff )#"Homo_sapiens.GRCh37.70.gtf" 
    halfwinwidth = int(sys.argv[3])
    #fragmentsize = 36

    ##check if a read overlaps a window
    #firepos = HTSeq.GenomicArrayOfSets( "auto", stranded=False )

    firepos = defaultdict(list)
    featurelist = []		        
    for feature in gfffile:
	chrom = feature.iv.chrom
        start = str(feature.iv.start)
        end = str(feature.iv.end)
	summit = abs(feature.iv.start-feature.iv.end)/2 + min(feature.iv.start, feature.iv.end)
	featurelist.append(chrom+'_'+start+'_'+end+'_'+str(summit))
	if summit < halfwinwidth:
            fire = HTSeq.GenomicInterval(feature.iv.chrom, 0, summit+halfwinwidth, '+')
	    firewindow = HTSeq.GenomicArray("auto", stranded = True)
	    #for i in xrange(0,len(range(0, summit+halfwinwidth, (summit + halfwinwidth)/10))-1):
		#firewindowiv = HTSeq.GenomicInterval(feature.iv.chrom,range(0, summit+halfwinwidth, (summit + halfwinwidth)/10)[i],range(0, summit+halfwinwidth, (summit + halfwinwidth)/10)[i+1],'.')
	        #firewindowset = HTSeq.GenomicArrayOfSets("auto",stranded=False)
		#firewindow [ firewindowiv ] += 0
	    firepos [ fire ].append(firewindow)
        else:
	    fire = HTSeq.GenomicInterval(feature.iv.chrom, summit-halfwinwidth, summit+halfwinwidth, '+')
	    firewindow = HTSeq.GenomicArray("auto", stranded = True)
	    #for j in xrange(0,len(range(summit-halfwinwidth, summit+halfwinwidth, halfwinwidth/5))-1):
		#firewindowiv = HTSeq.GenomicInterval(feature.iv.chrom,range(summit-halfwinwidth, summit+halfwinwidth, halfwinwidth/5)[j],range(summit-halfwinwidth, summit+halfwinwidth, halfwinwidth/5)[j+1],'.')      
                #firewindow [ firewindowiv ] += summit-halfwinwidth
	    firepos [ fire ].append(firewindow)

    firepos = count(firepos,bamfile,halfwinwidth)
    
    mylist = []
    mylist = range(0,100)
    writer.writerows([mylist])

    for f in featurelist:
	fchrom = f.split('_')[0]
	fstart = int(f.split('_')[1])
	fend = int(f.split('_')[2]) 
	fsummit = int(f.split('_')[3])#abs(fstart-fend)/2 + min(fstart, fend)
        line = []
	line.append(f)
        for fire in firepos:
	    firesummit = abs(fire.start-fire.end)/2+min(fire.start,fire.end)#-halfwinwidth
	    if fchrom == fire.chrom and fsummit == firesummit:
	        for firewindow in firepos[fire]:
                    firewindownew = defFeature(firewindow,f,halfwinwidth)
		    if firesummit < halfwinwidth:
		        for i in xrange(0,len(range(0, firesummit+halfwinwidth, (firesummit + halfwinwidth)/100))):
                            firewindowiv = HTSeq.GenomicInterval(fire.chrom,range(0, firesummit+halfwinwidth, (firesummit + halfwinwidth)/100)[i],range(0, firesummit+halfwinwidth, (firesummit + halfwinwidth)/100)[i]+(firesummit + halfwinwidth)/100,'-') 
	                    for iv, val in firewindownew[firewindowiv].steps():
	                        #line.append(val/42.187318)#chipWCE
				#line.append(val/29.748215)#chipNanog
				#line.append(val/31.706034)#chipSox2
				#line.append(val/28.759245)#chipOct4
				#line.append(val/1.829480)#groseqRep1
				#line.append(val/0.604995)#groseqRep2
				line.append(val/0.715728)#groseqRep3
		    else:
			for j in xrange(0,len(range(firesummit-halfwinwidth, firesummit+halfwinwidth, halfwinwidth/50))):
                            firewindowiv = HTSeq.GenomicInterval(fire.chrom,range(firesummit-halfwinwidth, firesummit+halfwinwidth, halfwinwidth/50)[j],range(firesummit-halfwinwidth, firesummit+halfwinwidth, halfwinwidth/50)[j]+halfwinwidth/50,'-')
			    for iv, val in firewindownew[firewindowiv].steps():
				#line.append(val/42.187318)#
				#line.append(val/29.748215)#
				#line.append(val/31.706034)#
				#line.append(val/28.759245)##normalize with total read counts
				#line.append(val/1.829480)#groseqRep1
                                #line.append(val/0.604995)#groseqRep2
                                line.append(val/0.715728)#groseqRep3
        #print len(line)
        writer.writerows([line])

    ofile.close()

if __name__=='__main__':
    sys.exit(main(sys.argv))

