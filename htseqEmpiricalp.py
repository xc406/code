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

def poisson_probability(nummotif, mlambda):
    # naive:   math.exp(-mlambda) * mlambda**nummotif / factorial(nummotif)
    # iterative, to keep the components from getting too large or small:
    p = math.exp(-mlambda)
    for i in xrange(nummotif):
        p *= mlambda
        p /= i+1
    return p

def main(argv):
    if len(argv) < 3:
        sys.stderr.write("Usage: %s htseq_output_file htseq_output_file_random\n" % argv[0])
        #sys.stderr.write("Usage: %s htseq_output_file\n" % argv[0])
	return 1
    if not os.path.isfile(argv[1]):
        sys.stderr.write('Error: htseq_output_file %r was not found!\n' % argv[1])
        return 1
    if not os.path.isfile(argv[2]):
        sys.stderr.write('Error: htseq_output_file_random %r was not found!\n' % argv[2])
        return 1
    #if not os.path.isfile(argv[3]):
        #sys.stderr.write('Error: htseq_output_fileR %r was not found!\n' % argv[3])
        #return 1
    
    infile = sys.argv[1]
    infiler = sys.argv[2]
    ifile = open(infile,'rt')
    ifiler = open(infiler,'rt')
    reader = csv.reader(ifile, delimiter = '\t') 
    readerr = csv.reader(ifiler,delimiter = '\t')
    (path,fname) = os.path.split(infile)

    ofile = open(os.path.join(path, fname + '_empiricalp'), 'w')
    writer = csv.writer(ofile, delimiter = '\t')

#sortedlist = sorted(reader, key=operator.itemgetter(3), reverse=False)	    

    countlist = []
    #ppoisson = defaultdict(float)    

    for row in readerr:
	if not row[0] in ['no_feature','ambiguous','too_low_aQual','not_aligned','alignment_not_unique']:
	    #site = row[0]
	    #countdict[site].append(int(row[1])) ## create a dictionary of gname and a list of read counts for all motif targets around each gene
	    countlist.append(int(row[1]))
    #print len(countlist)
    #print countdict
    for row in reader:	
	if not row[0] in ['no_feature','ambiguous','too_low_aQual','not_aligned','alignment_not_unique']:
	    #print sum(1 for c in countlist if c > int(row[1]))+1
            row.append(-1*math.log10((sum(1 for c in countlist if c > int(row[1]))+1.0)/(len(countlist)+1.0)))
	    #print row
            writer.writerows([row])
    
    ofile.close()

if __name__=='__main__':
    sys.exit(main(sys.argv))



