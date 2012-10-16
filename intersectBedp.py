import sys
import os
import csv
#import operator

def main(argv):
    if len(argv) < 2:
        sys.stderr.write("Usage: %s intersect_bed_file\n" % argv[0])
        return 1
    if not os.path.isfile(argv[1]):
        sys.stderr.write('Error: intersect_bed_file %r was not found!\n' % argv[1])
        return 1

    #ofile = open('/home/xc406/data/bed_test/' +  shortname + '_o', 'w')
    #writer = csv.writer(ofile, delimiter = '\t')    
    infile = sys.argv[1]

    (path,fname) = os.path.split(infile)
    (shortname, extension) = os.path.splitext(fname)
    ifile = open(infile,'rU')
    reader = csv.reader(ifile, delimiter = '\t')
    
    ofile = open('/home/xc406/data/bed_test/' +  shortname + '_o', 'w')
    writer = csv.writer(ofile, delimiter = '\t')

    i = 0
    gname = ''
    for row in reader:
	if i % 2 == 0:
	    line = row[-1].split(' ')
	    gname = line[1][:-1]
	else:
	    newline = ['', 0]
	    newline[1] = row[-3]
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

if __name__=='__main__':
    sys.exit(main(sys.argv))



