import sys
import os
import csv
import random

def main(argv):
    if len(argv) < 2:
        sys.stderr.write("Usage: %s upstream10kb_bed_file\n" % argv[0])
        return 1
    if not os.path.isfile(argv[1]):
        sys.stderr.write('Error: upstream10kb_bed_file %r was not found!\n' % argv[1])
        return 1
##    if not os.path.isfile(argv[2]):
##        sys.stderr.write('Error: TF_Info_file %r was not found!\n' % argv[2])
##        return 1
##    if not os.path.isfile(argv[3]):
##        sys.stderr.write('Error: nm#_conversion_file %r was not found!\n' % argv[3])
##        return 1
    
    infile = sys.argv[1]

    (path,fname) = os.path.split(infile)

    ifile = open(infile,'rt')
    reader = csv.reader(ifile, delimiter = '\t')
    ##writer = csv.writer(ofile, delimiter = '\t')

    ##(path1,fname1) = os.path.split(infile1)
    ##(path2,fname2) = os.path.split(infile2)

    #ifile = open('/home/xc406/data/mm9upstream10kb.bed','rt')
    #reader = csv.reader(ifile, delimiter = '\t')

    ofile = open('/home/xc406/data/background/Stat3_random_new.gff', 'w')
    writer = csv.writer(ofile, delimiter = '\t')
    
    bedlist = []
    for row in reader:
	if not '_' in row[0]:
	    bedlist.append((row[0],int(row[1]),row[-1],row[3]))
    	    ## a list consists of [chr, start, strand, nmnumber]
	
    #bedlist = list([(row[0],int(row[1]),row[-1],row[3]) for row in reader])
    #print bedlist
    #locilist = []
    for i in range(186794):
	loci = random.choice(bedlist)
	#print loci[0]
	pos = random.choice(range(loci[1],(loci[1]+9792)))
	#print pos
    	start = str(pos)
	stop = str(pos+208)
	row = [0,0,0,0,0,0,0,0,0]
	row[0] = loci[0]#chr
	line = loci[-1].split('_')
	row[1] = line[0] + '_' + line[1]#transcript_id
	row[2] = 'random'#feature
	row[3] = start#start
	row[4] = stop#stop
	row[5] = 1#pval
	row[6] = loci[2]#strand
	row[7] = '.'#frame
	row[8] = 'gene_id ' + row[1] + '; ' 'transcript_id ' + row[1]#group 
	writer.writerows([row])
    
    ofile.close()


if __name__=='__main__':
    sys.exit(main(sys.argv))




