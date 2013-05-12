### this script format htseq_output for aupr calculation
import sys
import os
import csv
#import operator
from collections import defaultdict

def main(argv):
    if len(argv) < 2:
        sys.stderr.write("Usage: %s htseq_output_file\n" % argv[0])
        return 1
    if not os.path.isfile(argv[1]):
        sys.stderr.write('Error: htseq_output_file %r was not found!\n' % argv[1])
        return 1

    #ofile = open('/home/xc406/data/bed_test/' +  shortname + '_o', 'w')
    #writer = csv.writer(ofile, delimiter = '\t')    
    infile = sys.argv[1]

    (path,fname) = os.path.split(infile)
    (shortname, extension) = os.path.splitext(fname)
    ifile = open(infile,'rU')
    reader = csv.reader(ifile, delimiter = '\t')
    
    ofile = open('/home/xc406/data/htseq_output/' +  shortname + '_gname.txt', 'wt')
    writer = csv.writer(ofile, delimiter = '\t')

    countdict = defaultdict(int)
    line = ['','',0,0]## store gs as [target_gene,chr,start,count]
    for row in reader:
	if not row[0] in ['no_feature','ambiguous','too_low_aQual','not_aligned','alignment_not_unique']:
	    countdict[row[0].split('_')[0]] += int(row[1]) 

    for gname in countdict:	    
        line = [gname,countdict[gname]]
	writer.writerows([line])    

    ofile.close()

if __name__=='__main__':
    sys.exit(main(sys.argv))



