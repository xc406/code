import sys
import os
import csv
from collections import defaultdict

def main(argv):
    if len(argv) < 3:
        sys.stderr.write("Usage: %s exp_fimo_output_file bkgd_fimo_output_file\n" % argv[0])
        return 1
    if not os.path.isfile(argv[1]):
        sys.stderr.write('Error: experiment_fimo_output_file %r was not found!\n' % argv[1])
        return 1
    if not os.path.isfile(argv[2]):
        sys.stderr.write('Error: background_fimo_output_file %r was not found!\n' % argv[2])
        return 1
##    if not os.path.isfile(argv[3]):
##        sys.stderr.write('Error: nm#_conversion_file %r was not found!\n' % argv[3])
##        return 1
    
    infile_exp = sys.argv[1]
    infile_bkgd = sys.argv[2]

    (path,fname) = os.path.split(infile_exp)
    (n,ext) = os.path.splitext(fname)
    ifile_exp = open(infile_exp,'rt')
    reader_exp = csv.reader(ifile_exp, delimiter = '\t')
    
    ifile_bkgd = open(infile_bkgd, 'rt')
    reader_bkgd = csv.reader(ifile_bkgd, delimiter = '\t')
    
    ofile = open('/home/xc406/data/mtl/' + n +'_count', 'wt')
    writer = csv.writer(ofile, delimiter = '\t')

    ##(path1,fname1) = os.path.split(infile1)
    ##(path2,fname2) = os.path.split(infile2)

    ##ifiletest = open('/home/xc406/data/mtl/cluster1fimoout_100812.txt','rt')
    ##readertest = csv.reader(ifiletf, delimiter = '\t')

    ## store motif counts per (random) target site in motifcount dictionaries 
    modict = defaultdict(list)
    #targetlist = []
    #tempmo = ''
    #temptgt = ''
    tgtmocount = defaultdict(int)
    for row in reader_bkgd:
	if '_0.62' in row[0]:
	    #print row
	    #if not (row[3],row[0]) in motifcount.keys():
            ## create a dict with keys as (target id ,motif id) and values as counts
	    tgtmocount[(row[1],row[0])] += 1 
    
    #print tgtmo	
    #ifile_bkgd.seek(0)
    #for tm in tgtmo:
	#for row in reader_bkgd:
	    #if row[0] == tm[0] and row[1] == tm[1]:		    
    	    	#tgtmocount[tm] += 1
	    #tgtmocount[]
	    
    #print tgtmocount
    while tgtmocount:
	
	(tgt, mo), cnt = tgtmocount.popitem()
	#print tgt, mo, cnt
	modict[mo].append(cnt)   
	
    print modict

    expmocount = defaultdict(int)
    for row in reader_exp:
	if '_0.62' in row[0]:
	    if row[1] == 'Tcra':
	    	expmocount[row[0]] += 1

    counter = defaultdict(int)
    pval = defaultdict(int)
    for mo in modict.keys():
	for cnt in modict[mo]:
	    if cnt >= expmocount[mo]:
		counter[mo] += 1
	## calculate empirical p-val
	pval[mo] = counter[mo]/(length(modict[mo])+1)
	
    for e in pval.keys():
	newrow = [e, pval[e]]
	writer.writerows([newrow])
    	
    ofile.close()
	
if __name__=='__main__':
    sys.exit(main(sys.argv))




