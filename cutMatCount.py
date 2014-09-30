import HTSeq
from matplotlib import pyplot
import os, sys, csv
from os import getenv
import pickle
from collections import defaultdict
import numpy as np

"""store featured reads in a [genomic feature x experiment] numpy array"""
def count(cutMat, featurelist, bamfile, ct):
    """count the number of reads within each bin/genomic region/feature"""
    i = 0.0
    for f in featurelist:
	chrombin = HTSeq.GenomicArray("auto", stranded = True)
        fchrom = f.split('_')[0]
        fstart = int(f.split('_')[1])
        fend = int(f.split('_')[2])
	fstrand = str(f.split('_')[3])
	fsize = abs(fend-fstart)
	#print f
	bwindowiv = HTSeq.GenomicInterval(fchrom, fstart, fend, fstrand)
        for almnt in bamfile[bwindowiv]:
	    #print almnt
	    #if not almnt is None:
            if almnt.aligned:
		#print almnt
		#if almnt.iv.strand == fstrand:
	        chrombin [ bwindowiv ] += 1
 	    #print almnt
        bwindowCount = 0.0
	for iv, val in chrombin[bwindowiv].steps():
	    #print iv, val, fsize
	    bwindowCount += float(val)#*200/fsize
	cutMat[i,ct] = bwindowCount
	i += 1 
    return cutMat

def main(argv):
    if len(argv) < 2:
        sys.stderr.write("Usage: %s \n"
				    "bam_file\n"
                                    "gff_file\n" % argv[0])
        return 1
    if not os.path.isfile(argv[1]):
        sys.stderr.write('Error: bam_file %r was not found!\n' % argv[1])
        return 1
    if not os.path.isfile(argv[2]):
        sys.stderr.write('Error: bed_file %r was not found!\n' % argv[8])
        return 1

    infile_bam = sys.argv[1]
    infile_gff = sys.argv[2]

    if getenv('MYTMP'):
	opath = getenv('MYTMP')
	print opath, 'output path'

    (bampath,bamfname) = os.path.split(infile_bam)
    (ctname,ctext) = os.path.splitext(bamfname)
    (gffpath,gfffname) = os.path.split(infile_gff)
    (tfname,tfext) = os.path.splitext(gfffname)

    #ofile = open(os.path.join(gffpath,'faireMat'),'wt')
    #writer = csv.writer(ofile, delimiter = '\t')

    bamfile1rep1 = HTSeq.BAM_Reader( infile_bam )

    #bamfile = HTSeq.BAM_Reader("wgEncodeUwDnaseTh17AlnRep1.bam") 

    gfffile = HTSeq.GFF_Reader( infile_gff )

    chrombin = HTSeq.GenomicArray("auto", stranded = False)
    featurelist = []		        
    for feature in gfffile:
	chrom = feature.iv.chrom
        start = str(feature.iv.start)#gff_reader shifted the base automatically for bam output 0-based
        end = str(feature.iv.end)
	strand = str(feature.iv.strand)
	#size = abs(feature.iv.end - feature.iv.start)
	#print start,end,size
	featurelist.append(chrom+'_'+start+'_'+end+'_'+strand)

    #faireMat = np.zeros((len(featurelist),10),dtype='float64')
    cutMat = np.memmap(os.path.join(opath, tfname + ctname + '.npy'), dtype = 'float64', mode = 'w+', shape = (len(featurelist),1))
    print 'processing', ctname
    cutMat = count(cutMat,featurelist,bamfile1rep1,0)
    #print 'processing th2'
    #cutMat = count(cutMat,featurelist,bamfile1rep2,fragmentsize, 1)
    #print 'processing treg'
    #cutMat = count(cutMat,featurelist,bamfile2rep1,fragmentsize, 2)
    #print 'processing th17'
    #cutMat = count(cutMat,featurelist,bamfile2rep2,fragmentsize, 3)
    #print 'processing th1Rep1'
    #cutMat = count(cutMat,featurelist,bamfile3rep1,fragmentsize, 4)
    #print 'processing th1Rep2'
    #cutMat = count(cutMat,featurelist,bamfile3rep2,fragmentsize, 5)
    #print 'processing th1Rep3'
    #cutMat = count(cutMat,featurelist,bamfile4rep1,fragmentsize, 6)

    np.savetxt(os.path.join(bampath, tfname + ctname + '.txt'), cutMat, delimiter = '\t')
    #mylist = []
    #mylist = range(0,10)
    #writer.writerows([mylist])

    #chromdict = defaultdict(tuple)
    #with open('chromdict.txt','w') as f:
        #pickle.dump(chromdict,f)
    #with open('chromdict.txt','r') as f:
        #chromdict = pickle.load(f)

if __name__=='__main__':
    sys.exit(main(sys.argv))

