
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

def main(argv):
    if len(argv) < 2:
        #sys.stderr.write("Usage: %s htseq_output_fileC htseq_output_fileL htseq_output_fileR\n" % argv[0])
        sys.stderr.write("Usage: %s htseq_output_file\n" % argv[0])
	return 1
    if not os.path.isfile(argv[1]):
        sys.stderr.write('Error: htseq_output_fileC %r was not found!\n' % argv[1])
        return 1
    #if not os.path.isfile(argv[2]):
        #sys.stderr.write('Error: htseq_output_fileL %r was not found!\n' % argv[2])
        #return 1
    #if not os.path.isfile(argv[3]):
        #sys.stderr.write('Error: htseq_output_fileR %r was not found!\n' % argv[3])
        #return 1
    
    infile = sys.argv[1]
    ifile = open(infile,'rt')
    reader = csv.reader(ifile, delimiter = '\t') 
    
    (path,fname) = os.path.split(infile)

    fname2 = ('_').join(fname.split('_')[:-1])
    ofile = open(os.path.join('/home/xc406/hg19priors/v073113/', fname2 + '_poissoncdf'), 'w')
    writer = csv.writer(ofile, delimiter = '\t')

    #c = Context(prec=500, rounding=ROUND_HALF_DOWN)
    #setcontext(c)

#sortedlist = sorted(reader, key=operator.itemgetter(3), reverse=False)	    

    countdict = defaultdict(list)
    countlist = []
    meandict = defaultdict(int)
    motifnumdict = defaultdict(int)
    lambdadict = defaultdict(int) 
    ppoisson = defaultdict(float)    

    for row in reader:
	if not row[0] in ['no_feature','ambiguous','too_low_aQual','not_aligned','alignment_not_unique']:
	    gname = row[0].split('_')[0]
	    countdict[gname].append(float(row[1]))#(float(row[2])) ## create a dictionary of gname and a list of read counts/empirical p vals for all motif targets around each gene
	    countlist.append(float(row[1]))#(float(row[2]))

    #print countdict

    for gname in countdict:
	meancount = sum(countdict[gname])/float(len(countdict[gname]))
	# meancount = np.mean(countdict[gname])
	meandict[gname] = meancount ## dictionary that stores mean read counts per gene
        motifnumdict[gname] = len(countdict[gname]) ## dictionary that stores number of motif hits per gene
    
    #print meandict
    #print motifnumdict

    for gname in meandict:
	with precision(5000):
	#lambdadict[gname] = sum(1 for c in countlist if c > meandict[gname])
	    lambdadict[gname] = BigFloat(10000.0)*sum(1 for c in countlist if c > meandict[gname])/BigFloat(2700000000.0) #1870000000.0 ## calculate lambda per gene
	#print gname

    #print lambdadict

    line = ['',0.0]#,0.0,0.0]
    for gname in lambdadict:
    	#ppoisson[gname]=(math.pow(lambdadict[gname],motifnumdict[gname])*math.exp(-lambdadict[gname]))/math.factorial(motifnumdict[gname]) ##this takes forever!
	#ppoisson[gname] = poisson_probability(motifnumdict[gname],lambdadict[gname])
	#ppoisson[gname] = -1*math.log10(max(1-poisson.cdf(motifnumdict[gname],lambdadict[gname]),1e-300))
	#ppoisson[gname] = -1*math.log10(1-poisson.cdf(motifnumdict[gname],lambdadict[gname]))
	#print ppoisson[gname]
	line[0],line[1]= gname,poisson_cdf(motifnumdict[gname],lambdadict[gname])#,line[2],line[3] = gname,motifnumdict[gname],float(lambdadict[gname]),poisson_cdf(motifnumdict[gname],lambdadict[gname])#ppoisson[gname]
	#p.append(line[1])
	#print line
	writer.writerows([line])
   #ppoisson_sorted = sorted(ppoisson.iteriterms(), key=operator.itemgetter(1), reverse=False)
    
    ofile.close()

if __name__=='__main__':
    sys.exit(main(sys.argv))



