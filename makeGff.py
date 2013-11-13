import sys
import os
import csv, pickle
from collections import defaultdict 

##assign windowsize to a smaller value, eg. 100 for 100bp sliding window with fixed winsize 1kb
def writeRegions(chrom,start,end,featureid,winSize,line,writer,offset):
    """output binned genomic regions in gff format"""
    for i in xrange(start+offset,end,winSize):
        line[0] = chrom #chr
        line[1] = '%08d' % featureid
        featureid += 1
        if end % winSize == 1 and i == range(start+offset,end,winSize)[len(range(start+offset,end,winSize))-1]:
            line[3] = i #range(chromdictClean[c][0],chromdictClean[c][1],1000)[i] #start
            line[4] = end#range(chromdictClean[c][0],chromdictClean[c][1],1000)[i] + 1000 #end
        else:
            line[3] = i
            line[4] = i+199
        line[-1] = 'gene_id ' + line[0] + '_' + str(line[3]) + '_' + str(line[4])
        writer.writerows([line])
    return featureid

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

    ofile = open('/home/xc406/data/mm9genome200bp.gff', 'w')
    writer = csv.writer(ofile, delimiter = '\t')
    #chr 14: 90,084,281-121,389,890
    

    ############only need to generate chromdict once and just pickle it next time
    #chromdict = defaultdict(list)
    #flag = True
    #pchrom = 'chr0'
    #pend = 0
    #chromdict[pchrom].append(0)
    #while flag:
	#row = ifile.readline()	
	#if len(row) > 1:
    	    #row = row.split('\n')[0].split('\t')
	    #chrom = row[0]
	    #if pchrom != chrom:
	        #chromdict[chrom].append(int(row[1]))
	        #chromdict[pchrom].append(pend)
	    #pchrom = chrom
            #pend = int(row[2])
	#else:
	    #flag = False
	    #chromdict[chrom].append(pend)

    #print chromdict
    #chromdictClean = defaultdict(list)
    #for c in chromdict:
	#if not 'random' in c:
    	    #if not c == 'chr0':
    		#chromdictClean[c] = chromdict[c]
    #with open('hg19chromdict.txt','w') as f:
	#pickle.dump(chromdictClean,f)
    with open('mm9chromdict.txt', 'r') as f:
	chromdictClean = pickle.load(f)

    line = ['','','feature','','','1','+','.','']
    featureid = 1
    print chromdictClean
    offset = 0
    winSize = 200
    for c in chromdictClean:
	#if c == 'chr14':##exclude anomaly regions
	    #writeRegions(c,chromdictClean[c][0],90084281,featureid,1000,line,writer)
	    #writeRegions(c,121389890,chromdictClean[c][1],featureid,1000,line,writer)
	#else:
	featureid = writeRegions(c,chromdictClean[c][0],chromdictClean[c][1],featureid,winSize,line,writer,offset)
	
    ofile.close()


if __name__=='__main__':
    sys.exit(main(sys.argv))

