import sys
import os
import csv
import random

def main(argv):
    if len(argv) < 2:
        sys.stderr.write("Usage: %s cluster_file\n" % argv[0])
        return 1
    if not os.path.isfile(argv[1]):
        sys.stderr.write('Error: cluster_file %r was not found!\n' % argv[1])
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

    ifile1 = open('/home/xc406/data/hg19.fa','rt')
    reader1 = csv.reader(ifile1, delimiter = '\t')

    #ofile = open('/home/xc406/data/mtl/cluster_background10.txt', 'w')
    #writer = csv.writer(ofile, delimiter = '\t')
    
    chromlist = []
    chromdict = {}
    for row in reader1:
	if '>' in row[0]:
            i = row[0][1:]
            if not i in chromlist:
	    	chromlist.append(i)
		chromdict[i] = 0
        if not '>' in row[0]:
	    chromdict[i] += len(row[0]) ## create a dictionary of chr and size
    		
    print chromdict
    poslist = []
    for row in reader:
	diff = int(row[2])-int(row[1])
	poslist.append(diff)
	#if not row[1] in chromlist:
	    #chromlist.append(row[1])
    
    #bedlist = list([(row[0],int(row[1]),row[-1],row[3]) for row in reader])
    #print bedlist
    #locilist = []
    for i in range(38000):
	chrom = random.choice(['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chrX','chrY'])
	loci = random.choice(poslist)
	pos1 = random.choice(range(1,(chromdict[chrom]-loci)))
	pos2 = pos1 + loci
	#newrow = ['', 0, 0, 0]
	newrow = [chrom, pos1, pos2, i+1]
	
	##print pos
    	#start = str(pos)
	#stop = str(pos+208)
	#row = [0,0,0,0,0,0,0,0,0]
	#row[0] = loci[0]#chr
	
	writer.writerows([newrow])
    
    ofile.close()

if __name__=='__main__':
    sys.exit(main(sys.argv))




