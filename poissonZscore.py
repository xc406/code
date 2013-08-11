
#calculate poisson give each gene target a poisson score based on empirical distribution of DHS read counts

##########
#implement scheme for calculating Ppoisson per motif window                
#                mean.pval.prox <- mean(M[ ix.prox ,7])
#                max.pval.prox <- max(M[ ix.prox ,7])
#                expected.num.peaks.genome.wide.prox <- length(which(M[,7]>=mean.pval.prox))
#                # lambda = num peaks / genome size * searched region
#                lambda.prox <- expected.num.peaks.genome.wide.prox/effective.genome.size*(tss.dist*2)
#                pval.pois.prox <- -log10(ppois(n.peaks.prox,lambda.prox,lower.tail=FALSE))
##########

import sys, os, csv
from collections import defaultdict
import math
import operator
#from scipy.stats import poisson
import numpy as np
#from decimal import *
from bigfloat import *

def poisson_cdf(nummotif, mlambda):
    # naive:   math.exp(-mlambda) * mlambda**nummotif / factorial(nummotif)
    # iterative, to keep the components from getting too large or small:
    with precision(5000):
        prob = BigFloat(1.0)
        for x in xrange(nummotif+1):
            p = exp(-mlambda)
	    for i in xrange(x+1):
	        if not i == 0:
                    p = mul(p,mlambda)
                    p = div(p, i)
	    prob = sub(prob, p)
	score = -log10(prob)

    return float(score)

def z_score(gname,plist,pdict):
    mu = np.mean(plist)
    stdev = np.std(plist)
    return (pdict[gname]-mu)/stdev

def main(argv):
    if len(argv) < 2:
        sys.stderr.write("Usage: %s poissonp_file\n" % argv[0])
	return 1
    if not os.path.isfile(argv[1]):
        sys.stderr.write('Error: poissonp_file %r was not found!\n' % argv[1])
        return 1
    
    infile = sys.argv[1]
    ifile = open(infile,'rt')
    reader = csv.reader(ifile, delimiter = '\t') 
    
    (path,fname) = os.path.split(infile)

    fname2 = ('_').join(fname.split('_')[:-1])
    ofile = open(os.path.join(path, fname + '_zscore'), 'w')
    writer = csv.writer(ofile, delimiter = '\t')

    pdict = defaultdict(float)
    plist = []
    #zdict = defaultdict(list)
    #zlist = []

    for row in reader:
	pdict[row[0]] = float(row[1]) 
	plist.append(float(row[1]))
    
    #for gname in pdict:
        #zdict[gname] = z_score(pdict[gname])

    #print zdict

    line = ['',0.0]#,0.0,0.0]
    for gname in pdict:
	line[0],line[1]= gname,z_score(gname,plist,pdict)#,line[2],line[3] = gname,motifnumdict[gname],float(lambdadict[gname]),poisson_cdf(motifnumdict[gname],lambdadict[gname])#ppoisson[gname]
	writer.writerows([line])

    ofile.close()

if __name__=='__main__':
    sys.exit(main(sys.argv))



