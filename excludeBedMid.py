import sys
import os
import csv

def main(argv):
    if len(argv) < 3:
        sys.stderr.write("Usage: %s bed_file TSS_bed_file \n" % argv[0])
        return 1
    if not os.path.isfile(argv[1]):
        sys.stderr.write('Error: bed_file %r was not found!\n' % argv[1])
        return 1
    if not os.path.isfile(argv[2]):
        sys.stderr.write('Error: TSS_bed_file %r was not found!\n' % argv[2])
        return 1

    infile = sys.argv[1]

    (path,fname) = os.path.split(infile)

    ifile = open(infile,'rt')
    ofile = open(os.path.join(path, "esc_nfr_ud2_mid.bed"),'wt')
    reader = csv.reader(ifile, delimiter = '\t')
    writer = csv.writer(ofile, delimiter = '\t')

    infile1 = sys.argv[2]

    ##(path1,fname1) = os.path.split(infile1)
    ##(path2,fname2) = os.path.split(infile2)

    ifile1 = open(infile1,'rt')
    reader1 = csv.reader(ifile1, delimiter = '\t')

    i = 0 
    #chrlist = []
    for row in reader:
	flag = 0
	ifile1.seek(0)
	for reg in reader1:
	    if reg[0] == row[0]:
		#if not row[0] in chrlist:
		    #chrlist.append(reg[0]) 
		#print 'chrom match', reg[0],row[0]
		#print range(int(reg[1]),int(reg[2])+1)
		#print (int(row[2])-int(row[1]))/2
		#print row[2],row[1]
		#if int(row[1]) > int(reg[1]) and int(row[1]) < int(reg[2])+1:
		    #print 'in, yes'
		mid = (int(row[2])-int(row[1]))/2 + int(row[1])
	        if mid in xrange(int(reg[1]),int(reg[2])+1): 
		#if int(row[1]) in xrange(int(reg[1]),int(reg[2])+1):
		    flag += 1
		    print 'mid-in', flag
		#elif int(row[2]) in xrange(int(reg[1]),int(reg[2])+1):
	            #flag += 1
		    #print 'right-in', flag
	#print row[1],row[2]
	if flag == 0:
	    writer.writerows([row])
	    i += 1
	    print i

    #print chrlist
    print 'final', i
    ofile.close()

if __name__=='__main__':
    sys.exit(main(sys.argv))
   
