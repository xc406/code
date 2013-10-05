import HTSeq, numpy
from matplotlib import pyplot
import os, sys, csv
import pickle
from collections import defaultdict
import numpy as np

"""store featured reads in a [genomic feature x experiment] numpy array"""

def count(faireMat, featurelist, bamfile, tot_reads, fragmentsize, ct):
    """count the number of reads within each bin/genomic region/feature"""
    i = 0.0
    #chrombin = HTSeq.GenomicArray("auto", stranded = False)
    for f in featurelist:
	chrombin = HTSeq.GenomicArray("auto", stranded = False)
        fchrom = f.split('_')[0]
        fstart = int(f.split('_')[1])
        fend = int(f.split('_')[2])
	if fstart == 0:
	    fstart = fragmentsize
	    print f
	bwindowiv = HTSeq.GenomicInterval(fchrom, fstart, fend, '.')
        for almnt in bamfile [bwindowiv]:
            if almnt.aligned:
		almnt.iv.length = fragmentsize
	        chrombin [ almnt.iv ] += 1
        bwindowCount = 0.0
	for iv, val in chrombin[bwindowiv].steps():
	    bwindowCount += val/tot_reads
	faireMat[i,ct] = bwindowCount
	i += 1 
    return faireMat

def main(argv):
    if len(argv) < 22:
        sys.stderr.write("Usage: %s bam_file1rep1\n"
				    "bam_file1rep2\n"
				    "bam_file2rep1\n"
				    "bam_file2rep2\n"
				    "bam_file3rep1\n"
				    "bam_file3rep2\n"
				    "bam_file4rep1\n"
				    "bam_file4rep2\n"
				    "bam_file_controlrep1\n"
				    "bam_file_controlrep2\n"	
                                    "gff_file\n"
				    "window_size\n"
				    "tot_reads1rep1\n"
				    "tot_reads1rep2\n"
				    "tot_reads2rep1\n"
				    "tot_reads2rep2\n"
				    "tot_reads3rep1\n"
				    "tot_reads3rep2\n"
				    "tot_reads4rep1\n"
				    "tot_reads4rep2\n"
				    "tot_reads_controlrep1\n"
				    "tot_reads_controlrep2\n" % argv[0])
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
        sys.stderr.write('Error: bam_file4rep2 %r was not found!\n' % argv[8])
        return 1
    if not os.path.isfile(argv[9]):
        sys.stderr.write('Error: bam_file_controlrep1 %r was not found!\n' % argv[9])
        return 1
    if not os.path.isfile(argv[10]):
        sys.stderr.write('Error: bam_file_controlrep2 %r was not found!\n' % argv[10])
        return 1    
    if not os.path.isfile(argv[11]):
        sys.stderr.write('Error: bed_file %r was not found!\n' % argv[11])
        return 1

    infile_bam1rep1 = sys.argv[1]
    infile_bam1rep2 = sys.argv[2]
    infile_bam2rep1 = sys.argv[3]
    infile_bam2rep2 = sys.argv[4]
    infile_bam3rep1 = sys.argv[5]
    infile_bam3rep2 = sys.argv[6]
    infile_bam4rep1 = sys.argv[7]
    infile_bam4rep2 = sys.argv[8]
    infile_bam_controlrep1 = sys.argv[9]
    infile_bam_controlrep2 = sys.argv[10]
    infile_gff = sys.argv[11]

    (bampath,bamfname) = os.path.split(infile_bam1rep1)
    
    (gffpath,gfffname) = os.path.split(infile_gff)
    (tfname,ext) = os.path.splitext(gfffname)

    #ofile = open(os.path.join(gffpath,'faireMat'),'wt')
    #writer = csv.writer(ofile, delimiter = '\t')

    bamfile1rep1 = HTSeq.BAM_Reader( infile_bam1rep1 )
    bamfile1rep2 = HTSeq.BAM_Reader( infile_bam1rep2 )
    bamfile2rep1 = HTSeq.BAM_Reader( infile_bam2rep1 )
    bamfile2rep2 = HTSeq.BAM_Reader( infile_bam2rep2 )
    bamfile3rep1 = HTSeq.BAM_Reader( infile_bam3rep1 )
    bamfile3rep2 = HTSeq.BAM_Reader( infile_bam3rep2 )
    bamfile4rep1 = HTSeq.BAM_Reader( infile_bam4rep1 )
    bamfile4rep2 = HTSeq.BAM_Reader( infile_bam4rep2 )
    bamfile_controlrep1 = HTSeq.BAM_Reader( infile_bam_controlrep1 )
    bamfile_controlrep2 = HTSeq.BAM_Reader( infile_bam_controlrep2 )
    #bamfile = HTSeq.BAM_Reader("wgEncodeUwDnaseTh17AlnRep1.bam") 

    gfffile = HTSeq.GFF_Reader( infile_gff )#"Homo_sapiens.GRCh37.70.gtf" 
    
    winwidth = int(sys.argv[12])

    tot_reads1rep1 = int(sys.argv[13])
    tot_reads1rep2 = int(sys.argv[14])
    tot_reads2rep1 = int(sys.argv[15])
    tot_reads2rep2 = int(sys.argv[16])
    tot_reads3rep1 = int(sys.argv[17])
    tot_reads3rep2 = int(sys.argv[18])
    tot_reads4rep1 = int(sys.argv[19])
    tot_reads4rep2 = int(sys.argv[20])
    tot_reads_controlrep1 = int(sys.argv[21])
    tot_reads_controlrep2 = int(sys.argv[22])
    fragmentsize = 36

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
    faireMat = np.memmap(os.path.join(gffpath,'faireMat.npy'), dtype = 'float64', mode = 'w+', shape = (len(featurelist),10))
    faireMat = count(faireMat,featurelist,bamfile1rep1,tot_reads1rep1,fragmentsize, 0)
    faireMat = count(faireMat,featurelist,bamfile1rep2,tot_reads1rep2,fragmentsize, 1)
    faireMat = count(faireMat,featurelist,bamfile2rep1,tot_reads2rep1,fragmentsize, 2)
    faireMat = count(faireMat,featurelist,bamfile2rep2,tot_reads2rep2,fragmentsize, 3)
    faireMat = count(faireMat,featurelist,bamfile3rep1,tot_reads3rep1,fragmentsize, 4)
    faireMat = count(faireMat,featurelist,bamfile3rep2,tot_reads3rep2,fragmentsize, 5)
    faireMat = count(faireMat,featurelist,bamfile4rep1,tot_reads4rep1,fragmentsize, 6)
    faireMat = count(faireMat,featurelist,bamfile4rep2,tot_reads4rep2,fragmentsize, 7)
    faireMat = count(faireMat,featurelist,bamfile_controlrep1,tot_reads_controlrep1,fragmentsize, 8)
    faireMat = count(faireMat,featurelist,bamfile_controlrep2,tot_reads_controlrep2,fragmentsize, 9)

    np.savetxt(os.path.join(gffpath, 'faireMat' + str(winwidth) + '.txt'), faireMat, delimiter = '\t')
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

