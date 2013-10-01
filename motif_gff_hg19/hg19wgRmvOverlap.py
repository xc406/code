import sys, os, csv
import bz2,time
#import operator

def modRow(prow, win):
    """fix motif size in the gff file and pad each motif window"""
    mlen = str(int(prow[-1].split(' ')[1].split('_')[3])+2)
    group = [prow[0],prow[3],prow[4],mlen,prow[6]]
    prow[-1] = 'gene_id ' + '_'.join(group) + '; sequence ' + prow[-1].split(' ')[-1]
    prow[3] = str(int(prow[3])-win)
    prow[4] = str(int(prow[4])+win)
    return prow

def minP(rowlist):
    """return the index of the row with the lowest pval"""
    pval = float(rowlist[0][5])
    minp = 0
    for p in xrange(len(rowlist)):
	if pval > float(rowlist[p][5]):
	    minp = p
	    pval = float(rowlist[p][5])
    return minp

def rmvOverlap(infile, ofile):
    #reader = mmap.mmap(ifile.fileno(),0,prot=mmap.PROT_READ)
    with bz2.BZ2File(infile,'r') as reader:
        #reader.seek(0)
        flag = True
	prow = []
	rowlist = []
	pstart,pend,psize,ppval = 0,0,0,1.0
	#i = 0
        while flag:
	    #i += 1
            row = reader.readline().split('\n')[0].split('\t')
            if not row is None:#i > 30:
                if len(row) > 1:
		    #print row, '<---original' 
		    start, end = int(row[3]), int(row[4])
		    size = end - start + 1
		    if pend in xrange(start,end+1):
               	        if abs(pend-start) > min(size, psize)/2.0:
                	    ##fimo_infile[motif_ID,chr,start,end,strand,score,pval,matched_seq]
                	    ##gff_outfile[chr,motif_ID,'motif',start,end,pval,strand,.,group(gene_id,gene_name,matched_seq)]
                            rowlist.append(prow)
			else:
			    if rowlist != []:
				rowlist.append(prow)
				#print rowlist
				wrow = rowlist[minP(rowlist)]
			        rowlist = []
				#print wrow, '<---best of the overlaps'
			    else:
				wrow = prow
			        #print wrow, '<---slightly overlap but all right'
                            ofile.writelines('\t'.join(modRow(wrow,100))+'\n')
		    else:
			if rowlist != []:
			    rowlist.append(prow)
			    #print rowlist
			    wrow = rowlist[minP(rowlist)]
			    rowlist = []
			    #print wrow, '<---best of the overlaps'
		            ofile.writelines('\t'.join(modRow(wrow,100))+'\n')
			elif len(prow) > 1:##just for the first entry to make sure [] does not get written
			    wrow = prow
			    #print wrow, '<---non-overlap'
		            ofile.writelines('\t'.join(modRow(wrow,100))+'\n')
		    prow = row
		    pstart, pend = start,end
                    psize = size
                else:
                    flag = False
            else:
                flag = False

def main(argv):
    if len(argv) < 2:
        sys.stderr.write("Usage: %s gff_file\n" % argv[0])
        return 1
    if not os.path.isfile(argv[1]):
        sys.stderr.write('Error: gff_file %r was not found!\n' % argv[1])
        return 1
	
    start_time = time.time()
    infile = sys.argv[1]

    (path,fname) = os.path.split(infile)
    #ofile = bz2.BZ2File('/home/xc406/data/final/' + fname, 'w')
    ofile = bz2.BZ2File('/scratch/xc406/hg19fimo/hg19gff1e3wg80/chrY/final/' + fname, 'w')

    rmvOverlap(infile,ofile)	

    ofile.close()
    print time.time() - start_time, "seconds"
	
if __name__=='__main__':
    sys.exit(main(sys.argv))



