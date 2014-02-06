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
#import pp
import numpy as np
import pickle
#from bigfloat import *

def main(argv):
    if len(argv) < 3:
        sys.stderr.write("Usage: %s emp_p_file htseq_output_file\n" % argv[0])
        #sys.stderr.write("Usage: %s htseq_output_file\n" % argv[0])
	return 1
    if not os.path.isfile(argv[1]):
        sys.stderr.write('Error: emp_p_file %r was not found!\n' % argv[1])
        return 1
    if not os.path.isfile(argv[2]):
        sys.stderr.write('Error: htseq_output_file %r was not found!\n' % argv[2])
        return 1
    
    infile = sys.argv[2]
    infiler = sys.argv[1]
    ifile = open(infile,'rt')
    ifiler = open(infiler,'rt')
    reader = csv.reader(ifile, delimiter = '\t') 
    readerr = csv.reader(ifiler,delimiter = '\t')
    (path,fname) = os.path.split(infile)

    ofile = open(os.path.join(path, fname+'Ll'), 'w')
    writer = csv.writer(ofile, delimiter = '\t')

#sortedlist = sorted(reader, key=operator.itemgetter(3), reverse=False)	    

    countdict = {}
    countlist = []
    #ppoisson = defaultdict(float)    

    for row in readerr:
	#if not row[0] in ['no_feature','ambiguous','too_low_aQual','not_aligned','alignment_not_unique']:
	count = float(row[0])
	if not row[1] == 'Inf':
	    countdict[count] = float(row[1])#,float(row[2])) ## create a dictionary of count and empirical pvalue
	    countlist.append(count)
    #print len(countlist)
    #print countdict
    for row in reader:
	#if not row[0] in ['no_feature','ambiguous','too_low_aQual','not_aligned','alignment_not_unique']:
	    #print sum(1 for c in countlist if c > int(row[1]))+1
	try:
	    row.append(countdict[round(float(row[-1]))])
	    #row.append(countdict[round(float(row[9]))][1])
	except KeyError:
	    #print "count not found", row[9]
	    row.append(countdict[max(countlist)])
	    #row.append(countdict[max(countlist)][1])
        writer.writerows([row])
    
    #with open('empVec.txt','w') as f:
        #pickle.dump(row,f)

    #with open('tsspos.txt','r') as f:
        #tsspos = pickle.load(f)

if __name__=='__main__':
    sys.exit(main(sys.argv))



