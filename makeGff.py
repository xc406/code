import sys
import os
import csv, pickle
from collections import defaultdict 

def writeRegions(chrom,start,end,featureid,windowsize,line,writer):
    """output binned genomic regions in gff format"""
    for i in xrange(start,end,windowsize):
        line[0] = chrom #chr
        line[1] = '%07d' % featureid
        featureid += 1
        if end % windowsize == 1 and i == range(start,end,windowsize)[len(range(start,end,windowsize))-1]:
            line[3] = i #range(chromdictClean[c][0],chromdictClean[c][1],1000)[i] #start
            line[4] = end#range(chromdictClean[c][0],chromdictClean[c][1],1000)[i] + 1000 #end
        else:
            line[3] = i
            line[4] = i+1000
        line[-1] = 'gene_id ' + line[0] + '_' + str(line[3]) + '_' + str(line[4])
        writer.writerows([line])
def main(argv):
    if len(argv) < 2:
        sys.stderr.write("Usage: %s bed_file\n" % argv[0])
        return 1
    if not os.path.isfile(argv[1]):
        sys.stderr.write('Error: bed_file %r was not found!\n' % argv[1])
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

    ofile = open('/home/xc406/data/esc/faire/genome1kb_chr14Test.gff', 'w')
    writer = csv.writer(ofile, delimiter = '\t')
    #chr 14: 90,084,281-121,389,890
    
#    chromdict = defaultdict(list)
#    flag = True
#    pchrom = 'chr0'
#    pend = 0
#    chromdict[pchrom].append(0)
#    while flag:
#	row = ifile.readline()	
#	if len(row) > 1:
#    	    row = row.split('\n')[0].split('\t')
#	    chrom = row[0]
#	    if pchrom != chrom:
#	        chromdict[chrom].append(int(row[1]))
#	        chromdict[pchrom].append(pend)
#	    pchrom = chrom
#            pend = int(row[2])
#	else:
#	    flag = False
#	    chromdict[chrom].append(pend)

    #print chromdict
#    chromdictClean = defaultdict(list)
#    for c in chromdict:
#	if not 'random' in c:
#    	    if not c == 'chr0':
#    		chromdictClean[c] = chromdict[c]
    #with open('chromdict.txt','w') as f:
	#pickle.dump(chromdictClean,f)
    with open('chromdict.txt', 'r') as f:
	chromdictClean = pickle.load(f)
    line = ['','','feature','','','1','+','.','']
    featureid = 1
    print chromdictClean
    for c in chromdictClean:
	if c == 'chr14':
	    writeRegions(c,chromdictClean[c][0],90084281,featureid,1000,line,writer)
	    #for i in xrange(chromdictClean[c][0],90084281,1000):
		#line[0] = c
		#line[1] = '%07d' % featureid
		#featureid += 1
		#if i == range(chromdictClean[c][0],90084281,1000)[len(range(chromdictClean[c][0],90084281,1000))-1]:
		    #line[3] = i
		    #line[4] = 90084281
	        #else:
		    #line[3] = i
		    #line[4] = i+1000
		#line[-1] = 'gene_id ' + line[0] + '_' + str(line[3]) + '_' + str(line[4])
		#writer.writerows([line])
	    writeRegions(c,121389890,chromdictClean[c][1],featureid,1000,line,writer)
	    #for i in xrange(121389890,chromdictClean[c][1],1000):
		#line[0] = c
		#line[1] = '%07d' % featureid
		#featureid += 1
	        #if chromdictClean[c][1] % 1000 == 1 and i == range(121389890,chromdictClean[c][1],1000)[len(range(121389890,chromdictClean[c][1],1000))-1]:
	            #line[3] = i
		    #line[4] = chromdictClean[c][1]
		#else:
		    #line[3] = i
                    #line[4] = i+1000
		#line[-1] = 'gene_id ' + line[0] + '_' + str(line[3]) + '_' + str(line[4])
		#writer.writerows([line])
	else:
	    writeRegions(c,chromdictClean[c][0],chromdictClean[c][1],featureid,1000,line,writer)
	    #for i in xrange(chromdictClean[c][0],chromdictClean[c][1],1000):
	        #line[0] = c #chr
	        #line[1] = '%07d' % featureid
	        #featureid += 1
	        #if chromdictClean[c][1] % 1000 == 1 and i == range(chromdictClean[c][0],chromdictClean[c][1],1000)[len(range(chromdictClean[c][0],chromdictClean[c][1],1000))-1]:
	            #line[3] = i #range(chromdictClean[c][0],chromdictClean[c][1],1000)[i] #start
	            #line[4] = chromdictClean[c][1]#range(chromdictClean[c][0],chromdictClean[c][1],1000)[i] + 1000 #end
	        #else:
		    #line[3] = i
		    #line[4] = i+1000
                #line[-1] = 'gene_id ' + line[0] + '_' + str(line[3]) + '_' + str(line[4])	 
                #writer.writerows([line])
    
    ofile.close()


if __name__=='__main__':
    sys.exit(main(sys.argv))




