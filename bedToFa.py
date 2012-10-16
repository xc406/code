import sys
import os
import csv
#import operator

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
    (shortname, extension) = os.path.splitext(fname)
    ifile = open(infile,'rt')
    reader = csv.reader(ifile, delimiter = '\t')
    ##writer = csv.writer(ofile, delimiter = '\t')

    ##(path1,fname1) = os.path.split(infile1)
    ##(path2,fname2) = os.path.split(infile2)

##ifile = open('/home/xc406/data/fimo052912.txt','rt')
##reader = csv.reader(ifile, delimiter = '\t')

    ofile = open('/home/xc406/data/' +  shortname + '_new.bed', 'w')
    writer = csv.writer(ofile, delimiter = '\t')
	
#sortedlist = sorted(reader, key=operator.itemgetter(3), reverse=False)	    

    ##remove all identicals and overlaps
    ifile.seek(0)
    #pline = ''
    pstart = ''
    pstrand = ''
    pchr = ''
    for row in reader:
	#if not row == pline:
	if not (row[0] == pchr and row[1] == pstart and row[-1] == pstrand):
	    pchr = row[0]
	    pstart = row[1]
	    pstrand = row[-1]
		
	    	#if not mstart in [int(pstart),int(pstop)]:
		    #if not  mstop in [int(pstart),int(pstop)]:
            writer.writerows([row])
	    #else:
		#writer.writerows([row])
	#pline = row
	#pstart = row[3]
        #pstop = row[4]
        #pchr = row[1]            

    ofile.close()

if __name__=='__main__':
    sys.exit(main(sys.argv))



