import sys
import os
import csv, pickle
from collections import defaultdict 

##process commandline output
##sort-bed
##bedops --chrom chr1 --closest --delim '\t' --dist sort-bed1.bed sort-bed2.bed > testALL 

ifile = open('/home/xc406/data/testALL','rt')
ofile = open('/home/xc406/data/chrALL','wt')

ifile.seek(0)

flag = True

while flag:
    line = ifile.readline()
    if line is None or len(line) < 2:
        flag = False
    else:
        seg = line.split('\r')
        r = seg[0].split('\t')[0:3]
        r.append(seg[1].split('\t')[-1])
        r.append(str(abs(int(seg[2].split('\t')[-1].split('\n')[0]))+4999)+'\n')
        row = '\t'.join(r)
        ofile.write(row)

ofile.close()

##assign windowsize to a smaller value, eg. 100 for 100bp sliding window with fixed winsize 1kb
def writeRegions(chrom,start,end,featureid,line,writer,gname,dist):
    """output binned genomic regions in gff format"""
    line[0] = chrom #chr
    line[1] = '%08d' % featureid
    featureid += 1
    #if end % winSize == 1 and i == range(start+offset,end,winSize)[len(range(start+offset,end,winSize))-1]:
    line[3] = str(int(start) + 1) ##gff 1-based! bed 0-based!#range(chromdictClean[c][0],chromdictClean[c][1],1000)[i] #start
    line[4] = str(int(end) + 1) #range(chromdictClean[c][0],chromdictClean[c][1],1000)[i] + 1000 #end
        #else:
            #line[3] = i
            #line[4] = i+199
    line[-1] = 'gene_id ' + str(gname) + '_' + line[0] + '_' + str(line[3]) + '_' + str(line[4]) + '_' + str(dist)
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

    ofile = open('/home/xc406/data/mm9genomegname.gff', 'w')
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
    #with open('mm9chromdict.txt', 'r') as f:
	#chromdictClean = pickle.load(f)

    line = ['','','feature','','','1','+','.','']
    featureid = 1
    #print chromdictClean
    offset = 0
    winSize = 200
    for row in reader:
	chrom,start,end,gname,dist = row[0],row[1],row[2],row[3],row[4]
	#if c == 'chr14':##exclude anomaly regions
	    #writeRegions(c,chromdictClean[c][0],90084281,featureid,1000,line,writer)
	    #writeRegions(c,121389890,chromdictClean[c][1],featureid,1000,line,writer)
	#else:
	featureid = writeRegions(chrom,start,end,featureid,line,writer,gname,dist)
	
    ofile.close()


if __name__=='__main__':
    sys.exit(main(sys.argv))

