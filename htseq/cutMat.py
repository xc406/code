import HTSeq, numpy
from matplotlib import pyplot
import os, sys, csv
import pickle
from collections import defaultdict
import numpy as np

"""store featured reads in a [genomic feature x experiment] numpy array"""
def count(cutMat, featurelist, bamfile, fragmentsize, ct):
    """count the number of reads within each bin/genomic region/feature"""
    i = 0.0
    #chrombin = HTSeq.GenomicArray("auto", stranded = False)
    for f in featurelist:
	chrombin = HTSeq.GenomicArray("auto", stranded = False)
        fchrom = f.split('_')[0]
        fstart = int(f.split('_')[1])
        fend = int(f.split('_')[2])
	#if fstart == 0:
	    #fstart = fragmentsize
	    #print f
	bwindowiv = HTSeq.GenomicInterval(fchrom, fstart, fend, '.')
        for almnt in bamfile: #[bwindowiv]:
	    #print almnt
            if almnt.aligned:
		#almnt.iv.length = fragmentsize
	        chrombin [ bwindowiv ] += 1
 	    #print almnt
        bwindowCount = 0.0
	for iv, val in chrombin[bwindowiv].steps():
	    bwindowCount += val
	cutMat[i,ct] = bwindowCount
	i += 1 
    return cutMat

def main(argv):
    if len(argv) < 16:
        sys.stderr.write("Usage: %s bam_file1rep1\n"
				    "bam_file1rep2\n"
				    "bam_file2rep1\n"
				    "bam_file2rep2\n"
				    "bam_file3rep1\n"
				    "bam_file3rep2\n"
				    "bam_file4rep1\n"
                                    "gff_file\n"
				    "window_size\n"
				    "tot_reads1rep1\n"
				    "tot_reads1rep2\n"
				    "tot_reads2rep1\n"
				    "tot_reads2rep2\n"
				    "tot_reads3rep1\n"
				    "tot_reads3rep2\n"
				    "tot_reads4rep1\n" % argv[0])
        return 1
    if not os.path.isfile(argv[1]):
        sys.stderr.write('Error: bam_file1rep1 %r was not found!\n' % argv[1])
        return 1
    if not os.path.isfile(argv[2]):
        sys.stderr.write('Error: bam_file1rep2 %r was not found!\n' % argv[2])
        return 1
    if not os.path.isfile(argv[3]):
        sys.stderr.write('Error: bam_file2rep1 %r was not found!\n' % argv[3])
        return 1
    if not os.path.isfile(argv[4]):
        sys.stderr.write('Error: bam_file2rep2 %r was not found!\n' % argv[4])
        return 1
    if not os.path.isfile(argv[5]):
        sys.stderr.write('Error: bam_file3rep1 %r was not found!\n' % argv[5])
        return 1
    if not os.path.isfile(argv[6]):
        sys.stderr.write('Error: bam_file3rep2 %r was not found!\n' % argv[6])
        return 1
    if not os.path.isfile(argv[7]):
        sys.stderr.write('Error: bam_file4rep1 %r was not found!\n' % argv[7])
        return 1
    if not os.path.isfile(argv[8]):
        sys.stderr.write('Error: bed_file %r was not found!\n' % argv[8])
        return 1

    infile_bam1rep1 = sys.argv[1]
    infile_bam1rep2 = sys.argv[2]
    infile_bam2rep1 = sys.argv[3]
    infile_bam2rep2 = sys.argv[4]
    infile_bam3rep1 = sys.argv[5]
    infile_bam3rep2 = sys.argv[6]
    infile_bam4rep1 = sys.argv[7]
    infile_gff = sys.argv[8]

    (bampath,bamfname) = os.path.split(infile_bam1rep1)
    
    (gffpath,gfffname) = os.path.split(infile_gff)
    (tfname,ext) = os.path.splitext(gfffname)

    #ofile = open(os.path.join(gffpath,'faireMat'),'wt')
    #writer = csv.writer(ofile, delimiter = '\t')

    bamfile1rep1 = HTSeq.SAM_Reader( infile_bam1rep1 )
    bamfile1rep2 = HTSeq.SAM_Reader( infile_bam1rep2 )
    bamfile2rep1 = HTSeq.SAM_Reader( infile_bam2rep1 )
    bamfile2rep2 = HTSeq.SAM_Reader( infile_bam2rep2 )
    bamfile3rep1 = HTSeq.SAM_Reader( infile_bam3rep1 )
    bamfile3rep2 = HTSeq.SAM_Reader( infile_bam3rep2 )
    bamfile4rep1 = HTSeq.SAM_Reader( infile_bam4rep1 )
    #bamfile = HTSeq.BAM_Reader("wgEncodeUwDnaseTh17AlnRep1.bam") 

    gfffile = HTSeq.GFF_Reader( infile_gff )#"Homo_sapiens.GRCh37.70.gtf" 
    
    winwidth = int(sys.argv[9])

    tot_reads1rep1 = int(sys.argv[10])
    tot_reads1rep2 = int(sys.argv[11])
    tot_reads2rep1 = int(sys.argv[12])
    tot_reads2rep2 = int(sys.argv[13])
    tot_reads3rep1 = int(sys.argv[14])
    tot_reads3rep2 = int(sys.argv[15])
    tot_reads4rep1 = int(sys.argv[16])
    fragmentsize = 1

    chrombin = HTSeq.GenomicArray("auto", stranded = False)
    featurelist = []		        
    for feature in gfffile:
	chrom = feature.iv.chrom
	#if not chrom == 'chrM':
        start = str(feature.iv.start+1)
        end = str(feature.iv.end)
	size = abs(feature.iv.end - feature.iv.start-1)
	featurelist.append(chrom+'_'+start+'_'+end+'_'+str(size))

    #faireMat = np.zeros((len(featurelist),10),dtype='float64')
    cutMat = np.memmap(os.path.join(gffpath,'cutMat.npy'), dtype = 'float64', mode = 'w+', shape = (len(featurelist),7))
    print 'processing cd4'
    cutMat = count(cutMat,featurelist,bamfile1rep1,fragmentsize, 0)
    print 'processing th2'
    cutMat = count(cutMat,featurelist,bamfile1rep2,fragmentsize, 1)
    print 'processing treg'
    cutMat = count(cutMat,featurelist,bamfile2rep1,fragmentsize, 2)
    print 'processing th17'
    cutMat = count(cutMat,featurelist,bamfile2rep2,fragmentsize, 3)
    print 'processing th1Rep1'
    cutMat = count(cutMat,featurelist,bamfile3rep1,fragmentsize, 4)
    print 'processing th1Rep2'
    cutMat = count(cutMat,featurelist,bamfile3rep2,fragmentsize, 5)
    print 'processing th1Rep3'
    cutMat = count(cutMat,featurelist,bamfile4rep1,fragmentsize, 6)

    np.savetxt(os.path.join(gffpath, 'cutMat' + str(winwidth) + '.txt'), faireMat, delimiter = '\t')
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

