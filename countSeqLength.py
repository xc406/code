###this script format luciferase assay gold standard
import sys
import os
import csv
#import operator
import pickle
from collections import defaultdict

def main(argv):
    if len(argv) < 2:
        sys.stderr.write("Usage: %s Fasta_file\n" % argv[0])
        return 1
    if not os.path.isfile(argv[1]):
        sys.stderr.write('Error: Fasta_file %r was not found!\n' % argv[1])
        return 1

    #ofile = open('/home/xc406/data/bed_test/' +  shortname + '_o', 'w')
    #writer = csv.writer(ofile, delimiter = '\t')    
    infile = sys.argv[1]

    (path,fname) = os.path.split(infile)
    (shortname, extension) = os.path.splitext(fname)
    ifile = open(infile,'rU')
    reader = csv.reader(ifile, delimiter = '\t')
    
    ofile = open('/home/xc406/data/esc/' +  shortname + '_size', 'wt')
    #writer = csv.writer(ofile, delimiter = '\t')

    #line = ['','',0,0,'',0.0,0.0]## store gs as [target_gene,chr,start,end,strand,ratio,fdr]
    #line = ['',0,0,'',0.0,'']## store gs as [chr,start,end,target_gene,ratio,strand]
    sizedict = defaultdict(int)
    size = 0
    flag = 0
    for row in reader:
	flag +=1
	#print row
	if not '>' in row[0]:
	    size += len(row[0])
	    #size = int(row[0].split('_')[2]) - int(row[0].split('_')[1]) 
	    #if size < 10:
	    #sizelist.append(size)
	else:
	    if flag != 1:
	        sizedict[seqname] = size
		ofile.write(str(size)+'\n')
	    seqname = row[0].split('>')[1]
	    size = 0
	#line = [row[5].split(':')[0],row[5].split(':')[1].split('-')[0],row[5].split(':')[1].split('-')[1],row[1],row[-2],row[-3]]

    with open('sizedict.txt','w') as f:
        pickle.dump(sizedict,f)

    print sizedict
    print 'Number of sequences ', len(sizelist)+1
    
    ofile.close()

if __name__=='__main__':
    sys.exit(main(sys.argv))



