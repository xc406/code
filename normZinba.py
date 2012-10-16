import sys
import os
import csv
import math
#import operator

def main(argv):
    if len(argv) < 2:
        sys.stderr.write("Usage: %s zinba_output_bed_file\n" % argv[0])
        return 1
    if not os.path.isfile(argv[1]):
        sys.stderr.write('Error: zinba_output_bed_file %r was not found!\n' % argv[1])
        return 1

    #ofile = open('/home/xc406/data/bed_test/' +  shortname + '_o', 'w')
    #writer = csv.writer(ofile, delimiter = '\t')    
    infile = sys.argv[1]

    (path,fname) = os.path.split(infile)
    (shortname, extension) = os.path.splitext(fname)
    ifile = open(infile,'rU')
    reader = csv.reader(ifile, delimiter = '\t')
    
    ofile = open('/home/xc406/data/FAIRE_seq/' +  shortname + '_logp.bed', 'w')
    writer = csv.writer(ofile, delimiter = '\t')

    for row in reader:
	newrow=['',0,0,'',0.0]
	#print row[4][:-1]
	newrow[0] = row[0]
	newrow[1] = row[1]
	newrow[2] = row[2]
	newrow[3] = 'zinba_peak'
	if float(row[4][:-1]) == 0.00000000000000:
	    newrow[4] = 310.0
	else:
	    newrow[4] = -1*math.log10(float(row[4][:-1]))
	 
	writer.writerows([newrow])
	#pline = row
	#pstart = row[3]
        #pstop = row[4]
        #pchr = row[1]            

    ofile.close()

if __name__=='__main__':
    sys.exit(main(sys.argv))



