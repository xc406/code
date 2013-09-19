import sys
import os
import csv
#from collections import defaultdict
import numpy as np
import bz2
import time
#import operator

start_time = time.time()

def map(list1,array):
    #a = np.array([])
    #size = int(list1[len(list1)-2].split('\t')[2])
    #a = np.zeros((size,))
    #pstart = int(list1[0].split('\t')[1])
    #a[0:pstart] = [0.0]*pstart
    #a = np.append(a,[0.0]*pstart)##set all non-recorded bps to 0    
    #print a
    liststrip = list1[:len(list1)-1]##get rid of the last empty entry
    for i in liststrip:
	#print i
        listtemp = []
        start = int(i.split('\t')[1])
        end = int(i.split('\t')[2])
        val = float(i.split('\t')[3])
        fragsize = end-start
        listtemp = [val]*fragsize
	array[start:end] = np.array(listtemp)
        #a = np.append(a,listtemp)#[val]*(end-start))
	#print array[start:end],np.array(listtemp)
    return array

def main(argv):
    if len(argv) < 3:
        sys.stderr.write("Usage: %s Mapability_file fimo_output_file (in bz2)\n" % argv[0])
        return 1
    if not os.path.isfile(argv[1]):
        sys.stderr.write('Error: Mapability_file %r was not found!\n' % argv[1])
	return 1
    if not os.path.isfile(argv[2]):
        sys.stderr.write('Error: fimo_output_file (in bz2) %r was not found!\n' % argv[2])
        return 1
    
    infile = sys.argv[1]
    infilefimo = sys.argv[2]

    (path,fname) = os.path.split(infilefimo)
    (shortname, ext) = os.path.splitext(fname)
    (mappath,mapfname) = os.path.split(infile)
    (mapshortname,mapext) = os.path.splitext(mapfname) 
    #ifile = open(infile,'rU')
    #ifilefimo = open(infilefimo,'rU')
 
    reader = bz2.BZ2File(infile, 'r')
    #readerfimo = bz2.BZ2File(infilefimo, 'r')
    ##writer = csv.writer(ofile, delimiter = '\t')

    ##(path1,fname1) = os.path.split(infile1)
    ##(path2,fname2) = os.path.split(infile2)

##ifile = open('/home/xc406/data/fimo052912.txt','rt')
##reader = csv.reader(ifile, delimiter = '\t')

    #ofile = open('/home/xc406/data/' + 'motif_mappability.gff', 'wt')
    #writer = csv.writer(ofile, delimiter = '\t')

    str1 = reader.read()
    list1 = str1.split('\n')
    mapchrom = list1[0].split('\t')[0]
    size = int(list1[len(list1)-2].split('\t')[2])
    #print size
    #a = np.zeros((size,), dtype = 'float32')
    #np.save(os.path.join(mappath,mapshortname), a)
    #print mapchrom
    ofile = bz2.BZ2File('/home/xc406/data/'+ shortname + 'map.bz2','w')
    #writer = csv.writer(ofile,delimiter = '\t')
    azeros = np.memmap(os.path.join(mappath,mapshortname+'.npy'), dtype = 'float64', mode = 'w+', shape = (size,))
    print len(azeros)
    amemmap = map(list1,azeros)
    #np.save(os.path.join(mappath,mapshortname), a)
    print len(amemmap)
    #amemmap = np.memmap(os.path.join(mappath,mapshortname+'.npy'), mode = 'r')
    with bz2.BZ2File(infilefimo,'r') as readerfimo:
    #for i in xrange(191417882):#len(readerfimo.readlines())): 
    	while True:
	    line = readerfimo.readline()
	    if not line is None:
	        row = line.strip().split('\t')
	        print row
	#chrom = row[1].split('_')[0]
	#if chrom == mapchrom:
	    #seqstart = int(row[1].split('_')[1])#for non-standard genome coordinates
	    #mstart,mend = seqstart+int(row[2])-1,seqstart+int(row[3])
                mstart,mend = int(row[2]),int(row[3])+1	
                print np.mean(amemmap[mstart:mend])
	        if np.mean(amemmap[mstart:mend]) > 0.8:
	    #print row, 'clean'
	            ofile.writelines(row)
	    else:
	        break

    #chromlist = []
    #chromdict = defaultdict(list)
    #for row in reader:
	#chrom = row[0]
	#start = int(row[1])
	#end = int(row[2])
        #chromdict[chrom][start:end-1] = [float(row[3])]*abs(start-end)

    #ifile.seek(0)
    print time.time() - start_time, "seconds"

    ofile.close()

if __name__=='__main__':
    sys.exit(main(sys.argv))



