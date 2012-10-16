import sys
import os
import csv
#import operator

def main(argv):
    if len(argv) < 4:
        sys.stderr.write("Usage: %s htseq_output_file1 htseq_output_file2 htseq_output_file3\n" % argv[0])
        return 1
    if not os.path.isfile(argv[1]):
        sys.stderr.write('Error: htseq_output_file1 %r was not found!\n' % argv[1])
        return 1
    if not os.path.isfile(argv[2]):
        sys.stderr.write('Error: htseq_output_file2 %r was not found!\n' % argv[2])
        return 1
    if not os.path.isfile(argv[3]):
        sys.stderr.write('Error: htseq_output_file3 %r was not found!\n' % argv[3])
        return 1
    
    infile = sys.argv[1]
    infile1 = sys.argv[2]
    infile2 = sys.argv[3]

    (path,fname) = os.path.split(infile)
    #(shortname, extension) = os.path.splitext(fname)
    ifile = open(infile,'rt')
    reader = csv.reader(ifile, delimiter = '\t')
    ifile1 = open(infile1,'rt')
    reader1 = csv.reader(ifile1, delimiter = '\t')
    ifile2 = open(infile2,'rt')
    reader2 = csv.reader(ifile2, delimiter = '\t')
    
##writer = csv.writer(ofile, delimiter = '\t')

    ##(path1,fname1) = os.path.split(infile1)
    ##(path2,fname2) = os.path.split(infile2)

##ifile = open('/home/xc406/data/fimo052912.txt','rt')
##reader = csv.reader(ifile, delimiter = '\t')

    ofile = open(path + fname + '_average', 'w')
    writer = csv.writer(ofile, delimiter = '\t')
	
#sortedlist = sorted(reader, key=operator.itemgetter(3), reverse=False)	    

    countdict = {}
    tot = 0.0
    tot1 =0.0
    tot2 = 0.0
    for row in reader:
	countdict[row[0]] = float(row[1])
        tot+=float(row[1])
    #glist = countdict.keys()
    #print countdict

    for k in countdict:
	countdict[k] /= tot

    for row in reader1:
	tot1+=float(row[1])
	#if not row[0] in glist:
	    #countdict[row[0]] = int(row[1])
	#else:
    ifile1.seek(0)
    for row in reader1:
	countdict[row[0]]+=float(row[1])/tot1
    
    #glist1 = countdict.keys()

    for row in reader2:
	tot2+=float(row[1])

    ifile2.seek(0)
    for row in reader2:
	#if not row[0] in glist:
        countdict[row[0]]+=float(row[1])/tot2
	new = ['',0]
	new[0] = row[0]
	new[1] = countdict[row[0]]*(tot+tot1+tot2)
	writer.writerows([new])
    print tot, tot1, tot2
    #print countdict
    ofile.close()

if __name__=='__main__':
    sys.exit(main(sys.argv))



