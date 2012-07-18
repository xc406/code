import sys
import os
import csv
#import operator

def main(argv):
    if len(argv) < 2:
        sys.stderr.write("Usage: %s gff_file\n" % argv[0])
        return 1
    if not os.path.isfile(argv[1]):
        sys.stderr.write('Error: gff_file %r was not found!\n' % argv[1])
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

##ifile = open('/home/xc406/data/fimo052912.txt','rt')
##reader = csv.reader(ifile, delimiter = '\t')

    ofile = open('/home/xc406/data/hg19gff_final/' +  fname, 'w')
    writer = csv.writer(ofile, delimiter = '\t')
	
#sortedlist = sorted(reader, key=operator.itemgetter(3), reverse=False)	    

    ##remove all identicals and overlaps
    ifile.seek(0)
    pline = ''
    pstart = ''
    pstop = ''
    pchr = 'chr1'
    for row in reader:
	if not row == pline:
	    if row[1] == pchr:
		mstart = int(row[3])
		mstop = int(row[4])
	    	if not ((mstart in [int(pstart),int(pstop)]) or (mstop in [int(pstart),int(pstop)])):
                    writer.writerows([row])
	    else:
		writer.writerows([row])
	pline = row
	pstart = row[3]
        pstop = row[4]
        pchr = row[1]            

    ofile.close()

if __name__=='__main__':
    sys.exit(main(sys.argv))



