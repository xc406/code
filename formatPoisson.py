###this script format Chip-Seq pPoisson and DGF pPoisson
import sys
import os
import csv
from collections import defaultdict

def main(argv):
    if len(argv) < 4:
        sys.stderr.write("Usage: %s chipPoisson_file dgfPoisson_file tss_file\n" % argv[0])
        return 1
    if not os.path.isfile(argv[1]):
        sys.stderr.write('Error: chipPoisson_file %r was not found!\n' % argv[1])
        return 1
    if not os.path.isfile(argv[2]):
        sys.stderr.write('Error: dgfPoisson_file %r was not found!\n' % argv[2])
        return 1
    if not os.path.isfile(argv[3]):
        sys.stderr.write('Error: tss_file %r was not found!\n' % argv[3])
        return 1

    #ofile = open('/home/xc406/data/bed_test/' +  shortname + '_o', 'w')
    #writer = csv.writer(ofile, delimiter = '\t')    
    chipfile = sys.argv[1]
    dgffile = sys.argv[2]
    tssfile = sys.argv[3]

    (path,fname) = os.path.split(chipfile)
    (shortname, extension) = os.path.splitext(fname)

    chipifile = open(chipfile,'rU')
    chipreader = csv.reader(chipifile, delimiter = '\t')
    
    dgfifile = open(dgffile,'rU')
    dgfreader = csv.reader(dgfifile, delimiter = '\t')

    tssifile = open(tssfile,'r')
    tssreader = csv.reader(tssifile, delimiter = '\t')

    ofile = open(os.path.join(path,shortname + '_llpoisson_se10_dgfcap'), 'wt')#'_gname_format'), 'wt')
    writer = csv.writer(ofile, delimiter = '\t')

    chippoissondict = defaultdict(float)
    dgfpoissondict = defaultdict(float)
    gnamelist = []

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
	dgfpoissondict[row[0]] = float(row[1])#(int(row[1]),float(row[2]))

    for row in tssreader:
	gnamelist.append(row[-1]) 

    print 'dgf->',len(dgfpoissondict),'chip->',len(chippoissondict),len(gnamelist)

    for gname in gnamelist:
	if (gname in chippoissondict) and (gname in dgfpoissondict):
            line = list([gname,chippoissondict[gname],dgfpoissondict[gname]])#[0],dgfpoissondict[gname][1]])
	    #line = list([gname,chippoissondict[gname],dgfpoissondict[gname][0],dgfpoissondict[gname][1]])
	elif (gname in chippoissondict) and (not gname in dgfpoissondict):
	    line = list([gname,chippoissondict[gname],0.0])
	elif (not gname in chippoissondict) and (gname in dgfpoissondict):
	    line = list([gname,0.0,dgfpoissondict[gname]])
	elif (not gname in chippoissondict) and (not gname in dgfpoissondict):
	    line = list([gname,0.0,0.0])
	    #line = list([gname,0.0,dgfpoissondict[gname][0],dgfpoissondict[gname][1]])
	writer.writerows([line])    

    #print max(p)
    #ofile.close()

if __name__=='__main__':
    sys.exit(main(sys.argv))



