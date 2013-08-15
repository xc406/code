###this script format Chip-Seq pPoisson and DGF pPoisson
import sys
import os
import csv
#import operator
from collections import defaultdict

def main(argv):
    if len(argv) < 3:
        sys.stderr.write("Usage: %s chipPoisson_file dgfPoisson_file\n" % argv[0])
        return 1
    if not os.path.isfile(argv[1]):
        sys.stderr.write('Error: chipPoisson_file %r was not found!\n' % argv[1])
        return 1
    if not os.path.isfile(argv[2]):
        sys.stderr.write('Error: dgfPoisson_file %r was not found!\n' % argv[1])
        return 1

    #ofile = open('/home/xc406/data/bed_test/' +  shortname + '_o', 'w')
    #writer = csv.writer(ofile, delimiter = '\t')    
    chipfile = sys.argv[1]
    dgffile = sys.argv[2]

    (path,fname) = os.path.split(chipfile)
    (shortname, extension) = os.path.splitext(fname)

    chipifile = open(chipfile,'rU')
    chipreader = csv.reader(chipifile, delimiter = '\t')
    
    dgfifile = open(dgffile,'rU')
    dgfreader = csv.reader(dgfifile, delimiter = '\t')

    ofile = open(os.path.join(path,shortname + '_emppoisson_dist_format'), 'wt')#'_gname_format'), 'wt')
    writer = csv.writer(ofile, delimiter = '\t')

    #p = []
    chippoissondict = defaultdict(float)
    dgfpoissondict = defaultdict(float)

    for row in chipreader:
        if not row[4] == 'strand':
	    #print row[9]
	    if not row[10] == 'NA': 
	        chippoissondict[row[0]] = float(row[10])
		#p.append(float(row[9]))
	    else:
		chippoissondict[row[0]] = 0.0
		
    for row in dgfreader:
	#dgfpoissondict[row[0]] = int(row[1])#
	dgfpoissondict[row[0]] = float(row[3])#(int(row[1]),float(row[2]))

    print len(dgfpoissondict),len(chippoissondict)

    for gname in dgfpoissondict:
	if gname in chippoissondict:
            line = list([gname,chippoissondict[gname],dgfpoissondict[gname]])#[0],dgfpoissondict[gname][1]])
	    #line = list([gname,chippoissondict[gname],dgfpoissondict[gname][0],dgfpoissondict[gname][1]])
	else:
	    line = list([gname,0.0,dgfpoissondict[gname]])#[0],dgfpoissondict[gname][1]])
	    #line = list([gname,0.0,dgfpoissondict[gname][0],dgfpoissondict[gname][1]])
	writer.writerows([line])    

    #print max(p)
    #ofile.close()

if __name__=='__main__':
    sys.exit(main(sys.argv))



