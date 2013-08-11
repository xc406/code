import sys
import os
import itertools
import csv
#import pickle
#import operator
from collections import defaultdict

def main(argv):
    if len(argv) < 3:
        sys.stderr.write("Usage: %s htseq_output_file %s centipede_gs\n" % argv[0])
        return 1
    if not os.path.isfile(argv[1]):
        sys.stderr.write('Error: htseq_output_file %r was not found!\n' % argv[1])
        return 1
    if not os.path.isfile(argv[2]):
        sys.stderr.write('Error: centipede_gs %r was not found!\n' % argv[2])
        return 1

    infile = sys.argv[1]
    gsfile=sys.argv[2]

    (path,fname) = os.path.split(infile)
    (gspath,gsfname) = os.path.split(gsfile)
    (gsname,gsext) = os.path.splitext(gsfname)

    ifile = open(infile,'rt')
    reader = csv.reader(ifile, delimiter = '\t')
 
    gsifile = open(gsfile,'rt')
    gsreader = csv.reader(gsifile,delimiter = '\t')

    hg19ud5file = open('/home/xc406/data/hg19ud5.bed','rt')
    hg19ud5reader = csv.reader(hg19ud5file,delimiter = '\t')

    ofile = open(os.path.join(gspath, gsname + '_format'), 'w')
    writer = csv.writer(ofile, delimiter = '\t')
	
#sortedlist = sorted(reader, key=operator.itemgetter(3), reverse=False)	    

    htseqdict = defaultdict(int)
    for row in reader:
	htseqdict[(row[1],int(row[2]))]= int(row[-1])

    #print htseqdict
    
    centdict = defaultdict(list)
    for row in gsreader:
	start = int(row[1])
	end = int(row[2])
        hg19ud5file.seek(0)
	for region in hg19ud5reader:
	    if region[0] == row[0]:##first check if gs is with 5kb tss
		if (start in range(int(region[1]),int(region[2]))) and (end in range(int(region[1]),int(region[2]))):  
	            centdict[(row[0],int(row[1]))] = [start, end]
		    break

    print 'number of gs', len(centdict)
	
    gsdict = defaultdict(list)
    for site in htseqdict:
	hit = 1
	for gs in centdict:
	    if site[0] == gs[0]:
		if site[1] in range(centdict[gs][0]-10,(centdict[gs][1]-centdict[gs][0])/2 + centdict[gs][0]):
		    gsdict[site] = [htseqdict[site],hit]
		    hit += 1
	
    print 'TP', len(gsdict)	

    for site in htseqdict:
	if not site in gsdict:
	    gsdict[site] = [htseqdict[site],0]

    print 'number of total predictions', len(gsdict), len(htseqdict)

    #with open('tgtdict.txt','w') as f:
        #pickle.dump(newtgtdict,f)	
        
    while gsdict:
	pred,match = gsdict.popitem()
	##gs_outfile[chr,pred_site,read_count,gs]
	newrow = list([pred[0],pred[1],match[0],match[1]]) 
	writer.writerows([newrow])
    ofile.close()
	

if __name__=='__main__':
    sys.exit(main(sys.argv))



