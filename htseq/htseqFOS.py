import sys
import os
import csv

def main(argv):
    if len(argv) < 4:
        sys.stderr.write("Usage: %s htseq_output_fileC htseq_output_fileL htseq_output_fileR\n" % argv[0])
        return 1
    if not os.path.isfile(argv[1]):
        sys.stderr.write('Error: htseq_output_fileC %r was not found!\n' % argv[1])
        return 1
    if not os.path.isfile(argv[2]):
        sys.stderr.write('Error: htseq_output_fileL %r was not found!\n' % argv[2])
        return 1
    if not os.path.isfile(argv[3]):
        sys.stderr.write('Error: htseq_output_fileR %r was not found!\n' % argv[3])
        return 1
    
    infileC = sys.argv[1]
    infileL = sys.argv[2]
    infileR = sys.argv[3]

    (path,fname) = os.path.split(infileC)
    #(shortname, extension) = os.path.splitext(fname)
    ifileC = open(infileC,'rt')
    readerC = csv.reader(ifileC, delimiter = '\t')
    ifileL = open(infileL,'rt')
    readerL = csv.reader(ifileL, delimiter = '\t')
    ifileR = open(infileR,'rt')
    readerR = csv.reader(ifileR, delimiter = '\t')
    
##writer = csv.writer(ofile, delimiter = '\t')

    ##(path1,fname1) = os.path.split(infile1)
    ##(path2,fname2) = os.path.split(infile2)

    ofile = open(os.path.join(path, fname + '_FOS'), 'w')
    writer = csv.writer(ofile, delimiter = '\t')
	
#sortedlist = sorted(reader, key=operator.itemgetter(3), reverse=False)	    

    countdict = {}
    countdictL = {}
    countdictR = {}
    totC = 0.0
    totL = 0.0
    totR = 0.0
    for row in readerC:
	countdict[row[0]] = float(row[1])+1
        totC+=float(row[1])
    #glist = countdict.keys()
    #print countdict

    #for k in countdict:
	#countdict[k] /= tot

    #for row in readerL:
	#totL+=float(row[1])
	#if not row[0] in glist:
	    #countdict[row[0]] = int(row[1])
	#else:
    ifileL.seek(0)
    for row in readerL:
	countdictL[row[0]] = countdict[row[0]]/(float(row[1])+1)
    
    #glist1 = countdict.keys()

    #for row in readerR:
	#totR+=float(row[1])

    ifileR.seek(0)
    for row in readerR:
	#if not row[0] in glist:
        countdictR[row[0]] =countdictL[row[0]] + countdict[row[0]]/(float(row[1])+1)
	new = ['',0]
	new[0] = row[0]
	new[1] = countdictR[row[0]]
	writer.writerows([new])
    
    #print countdict
    ofile.close()

if __name__=='__main__':
    sys.exit(main(sys.argv))



