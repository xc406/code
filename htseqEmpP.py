##########
#implement scheme for calculating Ppoisson per motif window                
#                mean.pval.prox <- mean(M[ ix.prox ,7])
#                max.pval.prox <- max(M[ ix.prox ,7])
#                expected.num.peaks.genome.wide.prox <- length(which(M[,7]>=mean.pval.prox))
#                # lambda = num peaks / genome size * searched region
#                lambda.prox <- expected.num.peaks.genome.wide.prox/effective.genome.size*(tss.dist*2)
#                pval.pois.prox <- -log10(ppois(n.peaks.prox,lambda.prox,lower.tail=FALSE))
##########

import sys
import os
import csv
from collections import defaultdict
import math
import operator
#from scipy.stats import poisson
import numpy as np
import pickle
from bigfloat import *

def main(argv):
    if len(argv) < 2:
        sys.stderr.write("Usage: %s htseq_output_file_random\n" % argv[0])
        #sys.stderr.write("Usage: %s htseq_output_file\n" % argv[0])
	return 1
    if not os.path.isfile(argv[1]):
        sys.stderr.write('Error: htseq_output_file_random %r was not found!\n' % argv[1])
        return 1
    #if not os.path.isfile(argv[3]):
        #sys.stderr.write('Error: htseq_output_fileR %r was not found!\n' % argv[3])
        #return 1
    
    #infile = sys.argv[1]
    infiler = sys.argv[1]
    #ifile = open(infile,'rt')
    ifiler = open(infiler,'rt')
    #reader = csv.reader(ifile, delimiter = '\t') 
    readerr = csv.reader(ifiler,delimiter = ' ')
    (path,fname) = os.path.split(infiler)

    ofile = open(os.path.join(path, 'empK562wg'), 'w')
    writer = csv.writer(ofile, delimiter = '\t')

#sortedlist = sorted(reader, key=operator.itemgetter(3), reverse=False)	    

    countlist = []
    #ppoisson = defaultdict(float)    

    for row in readerr:
	#if not row[0] in ['no_feature','ambiguous','too_low_aQual','not_aligned','alignment_not_unique']:
	    #site = row[0]
	    #countdict[site].append(int(row[1])) ## create a dictionary of gname and a list of read counts for all motif targets around each gene
	countlist.append(float(row[0]))
    #print len(countlist)
    #print countdict
    #for row in reader:	
    for num in xrange(1000):
	#if not row[0] in ['no_feature','ambiguous','too_low_aQual','not_aligned','alignment_not_unique']:
	    #print sum(1 for c in countlist if c > int(row[1]))+1
	row = [num]
	with precision(500):
            empP = -log10(div((sum(1 for c in countlist if c > int(num))+1.0),(len(countlist)+1.0)))
	    row.append(float(empP))
	print row
        writer.writerows([row])
    
    #with open('empVec.txt','w') as f:
        #pickle.dump(row,f)

    #with open('tsspos.txt','r') as f:
        #tsspos = pickle.load(f)

if __name__=='__main__':
    sys.exit(main(sys.argv))



