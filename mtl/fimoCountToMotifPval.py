import sys
import os
import csv
from collections import defaultdict
import pickle

def tgtmocnt(infile, version, tgtmodict):
    ifile = open(infile, 'rt')
    reader = csv.reader(ifile, delimiter = '\t')
    ifile.seek(0)
    for row in reader:
	if version in row[0]:
	    ## create a dictionary with (target id, motif id) as keys and motif counts as values 
	    tgtmodict[(row[1],row[0])] += 1
	#else:
	    #print 'wrong version of motif database :/'
    return tgtmodict

def mocnt(tgtmodict, modict):
    while tgtmodict:
	## create a dictionary of motif counts per target
	(tgt, mo), cnt = tgtmodict.popitem()
	modict[mo].append(cnt)
    return modict	
    
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
        
    ofile = open('/home/xc406/data/mtl/' + n +'tcrb_count', 'wt')
    writer = csv.writer(ofile, delimiter = '\t')

    ##(path1,fname1) = os.path.split(infile1)
    ##(path2,fname2) = os.path.split(infile2)
 
    modict = defaultdict(list)
    tgtmodict = defaultdict(int)
    #tgtmocountbkgd = tgtmocnt(infile_bkgd, '_0.62', tgtmodict)
    #with open('tgtmocountbkgd.txt','w') as f:
	#pickle.dump(tgtmocountbkgd,f)
    with open('tgtmocountbkgd.txt','r') as f:
	tgtmocountbkgd = pickle.load(f)
    #tgtmocountexp = tgtmocnt(infile_exp, '_0.62', tgtmodict)
    #with open('tgtmocountexp.txt','w') as f:
        #pickle.dump(tgtmocountexp,f)
    with open('tgtmocountexp.txt','r') as f:
        tgtmocountexp = pickle.load(f)
    #print tgtmocountbkgd

    #modictbkgd = mocnt(tgtmocountbkgd, modict)   
    #with open('modictbkgd.txt','w') as f:
        #pickle.dump(modictbkgd,f)
    with open('modictbkgd.txt','r') as f:
        modictbkgd = pickle.load(f)
    #print modictbkgd
    
    #modictexp = mocnt(tgtmocountexp, modict)
    #for mo in modictbkgd.keys():
        #print len(modictbkgd[mo])

    expmocount = defaultdict(int)
    ifile_exp.seek(0)
    for row in reader_exp:
	#if '_0.62' in row[0]:
	if row[1] == 'Tcrb':
	    expmocount[row[0]] += 1
    #print expmocount

    counter = defaultdict(int)
    pval = defaultdict(int)
    for mo in modictbkgd.keys():
	if mo in expmocount.keys():
	    #print expmocount[mo]
	    for cnt in modictbkgd[mo]:
	    	if cnt >= expmocount[mo]: #modictexp[mo]:
		    counter[mo] += 1
	else:
	    counter[mo] = 0
	print expmocount[mo]
	print modictbkgd[mo]
	print counter[mo]
	## calculate empirical p-val
	pval[mo] = (counter[mo]+1.0)/(len(modictbkgd[mo])+1)
    	
    for e in pval.keys():
	newrow = [e, pval[e]]
	writer.writerows([newrow])
    	
    ofile.close()
	
if __name__=='__main__':
    sys.exit(main(sys.argv))




