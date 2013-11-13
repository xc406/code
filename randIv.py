import sys, os, csv, random, pickle

"""generate random 10mers from genome and store 10mer genomic intervals in bed for BioStrings"""

def random_10mer(chromdict):
    """generate random chromosomal intervals size between 7 and 20"""
    chromlist = chromdict.keys()
    chrom = random.choice(chromlist)
    l = 10#random.randint(7, 20)
    start = random.randint(chromdict[chrom][0], chromdict[chrom][1])
    return (chrom, start, start+l)

def write_gff(nsample,chromdict,gfffile):
    gffwriter = csv.writer(gfffile, delimiter = '\t')
    for i in xrange(nsample):
        #loci = random.choice(bedlist)
        #print loci[0]
        #pos = random.choice(range(loci[1],(loci[1]+9790)))
        #print pos
        #start = str(pos)
        #stop = str(pos+212)
        chrom, start, end = random_10mer(chromdict)
        row = [0] * 9
        row[0] = chrom#chr
        #line = loci[-1].split('_')
        row[1] = 'M0000_0.80'#gene_id
        row[2] = 'random'#feature
        row[3] = start#start
        row[4] = end#end
        row[5] = 1#pval
        row[6] = random.choice(['+','-'])#strand
        row[7] = '.'#frame
        row[8] = 'gene_id ' + '%04d' % i + '_' + chrom + '_' + str(start) + '_' + str(end) + '_' + row[6] + '; seq'#group
        gffwriter.writerows([row])
    gfffile.close()

def write_bed(nsample,chromdict,bedfile):
    bedwriter = csv.writer(bedfile, delimiter = '\t')
    for i in xrange(nsample):
	chrom, start, end = random_10mer(chromdict)
	#row = [0] * 4
	row = chrom, start, end, '%04d' % i
	bedwriter.writerows([row])
    bedfile.close()

def main(argv):
    if len(argv) < 2:
        sys.stderr.write("Usage: %s number_of_samples\n" % argv[0])
        return 1
    #if not os.path.isfile(argv[1]):
        #sys.stderr.write('Error: ud5kb_bed_file %r was not found!\n' % argv[1])
        #return 1
##    if not os.path.isfile(argv[2]):
##        sys.stderr.write('Error: TF_Info_file %r was not found!\n' % argv[2])
##        return 1
##    if not os.path.isfile(argv[3]):
##        sys.stderr.write('Error: nm#_conversion_file %r was not found!\n' % argv[3])
##        return 1
    
    #infile = sys.argv[1]

    #(path,fname) = os.path.split(infile)

    #ifile = open(infile,'rt')
    #reader = csv.reader(ifile, delimiter = '\t')
    ##writer = csv.writer(ofile, delimiter = '\t')

    with open('hg19chromdict.txt', 'r') as f:
        chromdict = pickle.load(f)
    ##(path1,fname1) = os.path.split(infile1)
    ##(path2,fname2) = os.path.split(infile2)

    #ifile = open('/home/xc406/data/mm9upstream10kb.bed','rt')
    #reader = csv.reader(ifile, delimiter = '\t')

    #ofile = open('/home/xc406/data/myc_random3.gff', 'w')
    gfffile = open('/home/xc406/data/hg19random.gff', 'w')

    bedfile = open('/home/xc406/data/hg19random.bed', 'w')
  
    #bedlist = []
    #for row in reader:
	#if not '_' in row[0]:
	    #bedlist.append((row[0],int(row[1]),str(row[-1])))
    	    ## a list consists of [chr, start, gene_id]
	
    #bedlist = list([(row[0],int(row[1]),row[-1],row[3]) for row in reader])

    nsample = int(sys.argv[1])
    write_gff(nsample,chromdict,gfffile)
    write_bed(nsample,chromdict,bedfile)

if __name__=='__main__':
    sys.exit(main(sys.argv))




