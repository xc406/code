import sys
import os
import csv
import math
#import operator

def main(argv):
    if len(argv) < 2:
        sys.stderr.write("Usage: %s gff_file\n" % argv[0])
        return 1
    if not os.path.isfile(argv[1]):
        sys.stderr.write('Error: gff_file %r was not found!\n' % argv[1])
        return 1

    #ofile = open('/home/xc406/data/bed_test/' +  shortname + '_o', 'w')
    #writer = csv.writer(ofile, delimiter = '\t')    
    infile = sys.argv[1]

    (path,fname) = os.path.split(infile)
    (shortname, extension) = os.path.splitext(fname)
    ifile = open(infile,'r')
    reader = csv.reader(ifile, delimiter = '\t')
    
    ofile = open('/home/xc406/data/STAT3_hg19mappable.bed', 'wt')
    writer = csv.writer(ofile, delimiter = '\t')

    #pvalist = []
    #normpval = []
    #for row in reader:
	#print row
	#type(row[4])
	#newrow = ['',0,0,'',0]
	#newrow[0] = row[0]
	#newrow[1] = row[3]
	#newrow[2] = row[4]
	#newrow[3] = row[2]
	#logp = -10*math.log10(float(row[5]))
	#pvalist.append(logp)
	
    #minpval = min(pvalist)
    #maxpval = max(pvalist)
    #print minpval,maxpval
    #l = maxpval-minpval
    #for i in pvalist:
	#ni = ((i - minpval)/l)*3100.0
	#normpval.append(ni)
    #print len(normpval)
    #i=0
    ifile.seek(0)    
    for row in reader:
	newrow = ['','','','']
        newrow[0] = row[0]
        newrow[1] = row[3]
        newrow[2] = row[4]
        newrow[3] = row[5]
	#newrow[4] = normpval[i]
	#newrow[4] = str(pvalist[i])
	#i+=1
	writer.writerows([newrow])
       
	    #else:
		#writer.writerows([row])
	#pline = row
	#pstart = row[3]
        #pstop = row[4]
        #pchr = row[1]            

    ofile.close()

if __name__=='__main__':
    sys.exit(main(sys.argv))



