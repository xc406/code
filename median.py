import sys
import csv
import os
#from pylab import *
#import matplotlib.pyplot


def median(l):
    l.sort()
    if len(l)%2 == 0:
	i = len(l)/2
	median = (l[i-1] + l[i])/2.0
    else:
	i = len(l)/2
	median = l[i]
    return median

def motiflength(motifID, fimoerr):
    ##generate a dict of motifID and length for normalization            
    #global motifdict
    motifdict={}
    i=0
    for row in fimoerr:
        if i % 2 == 0:
            line = row[0].split(' ')
            motifID = line[2].lstrip('-')
            if not motifID in motifdict:
                motifdict[motifID] = line[-1][:-1]
        i+=1
    return motifdict[motifID]

    #motiflist=[]
    #motiflist = list([(row[0],row[-1]) for row in readerfile1])
    #motifdict = {}
    #for k,v in motiflist:
    #    if not k in motifdict:
    #        motifdict[k] = len(v)

    #print motifdict	

def read_count(htseq):
    tot_count = 0 #total number of reads in the DHS file
    for row in htseq:
        read = int(row[1])  
        tot_count += read
    return tot_count

def main(argv):
    if len(argv) < 4:
        sys.stderr.write("Usage: %s htseq-count_output_file gff_file fimo_stderr_file\n" % argv[0])
        return 1
    if not os.path.isfile(argv[1]):
        sys.stderr.write('Error: htseq-count_output_file %r was not found!\n' % argv[1])
        return 1
    if not os.path.isfile(argv[2]):
        sys.stderr.write('Error: gff_file %r was not found!\n' % argv[2])
        return 1
    if not os.path.isfile(argv[3]):
        sys.stderr.write('Error: fimo_stderr_file was not found!\n' % argv[3])
        return 1

    infile = sys.argv[1]
    infile1 = sys.argv[2]
    infile2 = sys.argv[3]
    
    (path,fname) = os.path.split(infile1)
    (shortname, extension) = os.path.splitext(fname)
    (patho,fname2) = os.path.split(infile)

    ifile = open(infile,'rt')
    ifile1 = open(infile1, 'rt')
    ifile2 = open(infile2, 'rt')

    #ofile = open(os.path.join(patho, shortname+ "_filter.gff"),'w')
    reader = csv.reader(ifile, delimiter = '\t')
    reader1 = csv.reader(ifile1, delimiter = '\t')
    reader2 = csv.reader(ifile2, delimiter = '\t')
    
##ofile = open('/Users/xichen/Downloads/Jund_filter.gff', 'w')
    #writer = csv.writer(ofile, delimiter = '\t')

##inputfile = open('/Users/xichen/Downloads/Jund_update_htseq_output1', 'r')
##readerfile = csv.reader(inputfile, delimiter = '\t')

##inputfile1 = open('/Users/xichen/Downloads/fimoout070312.err', 'rt')
##readerfile1 = csv.reader(inputfile1, delimiter = '\t')

##inputfile2 = open('/Users/xichen/Downloads/Jund_update.gff', 'r')
##readerfile2 = csv.reader(inputfile2, delimiter = '\t')

    nmlist = [] #a list of non-zero TF_names
    nmlistnormalized = []
    motifdict={}
    lenlist = []
    #previous = 8
    for row in reader2:
	line = row[0].split(' ')
        motifID = line[2].lstrip('-')
        if not motifID in motifdict:
            motifdict[motifID] = line[-1][:-1]
	    if not motifdict[motifID] == 'q-value':
	    	lenlist.append(int(motifdict[motifID]))

    print median(lenlist)	          
    #print motifdict
    #print lenlist
    ifile1.seek(0)
    for row in reader1:
	m = row[1]
	break
 
    ml = int(motiflength(m, reader2)) + 200
    ifile.seek(0)
    a = read_count(reader)
    #print a
    #print ml
   
    #ifile.seek(0)
    #ifile1.seek(0)
    #for row in reader:
	#print row
        #nmlistnormalized.append((float(row[1])/a/int(ml)))
    
    ifile.seek(0)
    for row in reader:
        read_num = float(row[1])
	read_normalized = read_num/a/ml
	#print read_normalized
	#if not read_normalized == 0.0:
        #if read_normalized > 3.0e-09:
	#if not read_num in [0.0, 1.0]:
	if read_normalized > 1.2e-10:
             nmlist.append(row[0])
    
    #print nmlist

    ifile1.seek(0)
    for row in reader1:
        line = row[-1].split(';')
        nm = line[0][8:-1]
        #print nm
        #nm = row[-3][:-1]
        if nm in nmlist:
             writer.writerows([row])
	     #print row
    
    ofile.close()

    #print tot_count
    #print len(nmlistnormalized)
    ###print len(nmlist)
    ####print nmlist
    ##figure()
    ##xVals = (nmlistnormalized)
    ##hist(xVals, bins = 100)
    ##show()

if __name__=='__main__':
    sys.exit(main(sys.argv))
