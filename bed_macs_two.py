#####################
##generate a three columns file with 
##gene names/ motif pval/ dhs macs pval 
##from intersected bed and macs output files
##and sort the output three columns file by gene names
#####################

import sys
import os
import csv
import math
import operator

def main(argv):
    if len(argv) < 2:
        sys.stderr.write("Usage: %s intersect_bed_file\n" % argv[0])
        return 1
    if not os.path.isfile(argv[1]):
        sys.stderr.write('Error: intersect_bed_file %r was not found!\n' % argv[1])
        return 1

    infile = sys.argv[1]

    (path,fname) = os.path.split(infile)
    (shortname, extension) = os.path.splitext(fname)
    ifile = open(infile,'rU')
    reader = csv.reader(ifile, delimiter = '\t')
    
    ofile = open('/home/xc406/data/mel_dhs/' +  shortname + '_two', 'w')
    writer = csv.writer(ofile, delimiter = '\t')

    i = 0
    gname = ''
    pval = ''
    for row in reader:
	if i % 2 == 0:
	    line = row[-1].split(' ')
	    gname = line[1][:-1]
	    pval = -1*math.log10(float(row[5]))
	else:
	    newline = ['', 0, 0]
	    newline[2] = float(row[-1])/10.0
	    newline[1] = pval
	    newline[0] = gname    
            writer.writerows([newline])
        i+=1

	    #else:
		#writer.writerows([row])
	#pline = row
	#pstart = row[3]
        #pstop = row[4]
        #pchr = row[1]            

    ofile.close()

    ifile2 = open('/home/xc406/data/mel_dhs/' +  shortname + '_two', 'r')
    reader2 = csv.reader(ifile2, delimiter = '\t')
    sortedlist = sorted(reader2, key = operator.itemgetter(0), reverse = False)
    
    ofile2 = open('/home/xc406/data/mel_dhs/' + shortname + '_two_psorted', 'w')
    writer2 = csv.writer(ofile2, delimiter = '\t')
    for item in sortedlist:
	writer2.writerows([item])
	
    ofile2.close()

if __name__=='__main__':
    sys.exit(main(sys.argv))



