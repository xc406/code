import sys
import os
import csv
from collections import defaultdict
import pickle
from math import log, exp, pow, floor

# global constants
_log10 = log(10)
_log_zero = -1e10               # Zero on the log scale.
_log_small = -0.5e10            # Threshold below which everything is zero.
_mm_nats = log(-_log_zero)      # dynamic range of double in NATs

def getLogFETPvalue(p, P, n, N, log_pthresh):
    """ Return log of hypergeometric pvalue of #pos >= p
            p = positive successes
            P = positives
            n = negative successes
            N = negatives
            log_pthresh = short-circuit if p will be greater
    """
    # check that p-value is less than 0.5
    # if p/float(P) > n/float(N):
    # if (p * N > n * P):
    if (p * N > n * P) and log_hyper_323(p,P,n+p,N+P) < log_pthresh:
        # apply Fisher Exact test (hypergeometric p-value)
        log_pvalue = log_getFETprob(N-n, n, P-p, p)[4];
    else:
        log_pvalue = 0          # pvalue = 1

    return log_pvalue


# Routines for computing the logarithm of a sum in log space.
def log_sum(logx, logy):
    """ Return the log(x+y) given log(x) and log(y). """
    if logx > logy:
        return log_sum1(logx, logy)
    else:
        return log_sum1(logy, logx)

def log_sum1(logx, logy):
    if (logx - logy) > _mm_nats:
        return logx
    else:
        return logx + log(1 + my_exp(logy - logx))

def my_exp(x):
    if x < _log_small:
        return 0.0
    else:
        return exp(x)

# print very large or small numbers
def sprint_logx(logx, prec, format):
    """ Print x with given format given logx.  Handles very large
    and small numbers with prec digits after the decimal.
    Returns the string to print."""
    log10x = logx/_log10
    e = floor(log10x)
    m = pow(10, (log10x - e))
    if ( m + (.5*pow(10,-prec)) >= 10):
        m = 1
        e += 1
    str = format % (m, e)
    return str

# Fisher's Exact Test
_log0_99999999 = log(0.9999999)
_log1_00000001 = log(1.00000001)
def log_getFETprob(a1, a2, b1, b2):
    """Computes Fisher's exact test based on a
    null-hypothesis distribution specified by the totals, and
    an observed distribution specified by b1 and b2, i.e.
    determines the probability of b's outcomes 1 and 2.

    Returns an immutable list consisting of the exact
    probability, and assorted p-values (sless, sright, sleft,
    slarg) based on the density."""
    (log_sless, log_sright, log_sleft, log_slarge) = (_log_zero, _log_zero, _log_zero, _log_zero)
    n = a1 + a2 + b1 + b2
    row1 = a1 + a2 # the row containing the null hypothesis
    col1 = a1 + b1 # the column containing samples for outcome 1
    max = row1
    if col1 < max:
        max = col1
    min = row1 + col1 - n
    if min < 0:
        min = 0
    if min == max:
        #rt = (prob, sless, sright, sleft, slarg) = (1.0,1.0,1.0,1.0,1.0)
        rt = (log_prob, log_sless, log_sright, log_sleft, log_slarg) = (0,0,0,0,0)
        return rt
    log_prob = log_hyper0(a1, row1, col1, n)
    log_sleft = _log_zero
    log_p = log_hyper(min)

    i = min + 1
    while log_p < (_log0_99999999 + log_prob):
        log_sleft = log_sum(log_sleft, log_p)
        log_p = log_hyper(i)
        i = i + 1

    i = i - 1
    if log_p < (_log1_00000001 + log_prob):
        log_sleft = log_sum(log_sleft, log_p)
    else:
        i = i - 1

    log_sright = _log_zero
    log_p = log_hyper(max)

    j = max - 1
    while log_p < (_log0_99999999 + log_prob):
        log_sright = log_sum(log_sright, log_p)
        log_p = log_hyper(j)
        j = j - 1

    j = j + 1
    if log_p < (_log1_00000001 + log_prob):
        log_sright = log_sum(log_sright, log_p)
    else:
        j = j + 1

    if abs(i - a1) < abs(j - a1):
        log_sless = log_sleft
        log_slarg = log_sum(1.0, -log_sleft)
        log_slarg = (log_slarg, log_prob)
    else:
        log_sless = log_sum(1.0, -log_sright)
        log_sless = (log_sless, log_prob)
        log_slarg = log_sright
    return (log_prob, log_sless, log_sright, log_sleft, log_slarg)

# log gamma function using continued fractions
def lngamm(z):
    x = 0.0
    x = x + 0.1659470187408462e-06/(z+7.0)
    x = x + 0.9934937113930748e-05/(z+6.0)
    x = x - 0.1385710331296526    /(z+5.0)
    x = x + 12.50734324009056     /(z+4.0)
    x = x - 176.6150291498386     /(z+3.0)
    x = x + 771.3234287757674     /(z+2.0)
    x = x - 1259.139216722289     /(z+1.0)
    x = x + 676.5203681218835     /(z)
    x = x + 0.9999999999995183
    return log(x)-5.58106146679532777-z+(z-0.5)*log(z+6.5)

# log n! computed using gamma function
_lnfact_hash = {}
def lnfact(n):
    if n<=1:
        return 0.0

    key = str(n)
    if _lnfact_hash.has_key(key):
        return _lnfact_hash[key]

    result = lngamm(n+1.0)
    _lnfact_hash[key] = result
    return result

# log binomial coefficient n choose k
def lnbico(n, k):
    return lnfact(n)-lnfact(k)-lnfact(n-k)

def log_hyper_323(n11, n1_, n_1, n):
    return lnbico(n1_,n11)+lnbico(n-n1_,n_1-n11)-lnbico(n,n_1)

_log_sprob = 0
def log_hyper0(n11i, n1_i, n_1i, ni):
    global _sn11, _sn1_, _sn_1, _sn, _log_sprob
    if not ((n1_i|n_1i|ni)!=0):
        if not (n11i % 10 == 0):
            if n11i==_sn11+1:
                _log_sprob = _log_sprob + \
                        log ( ((_sn1_-_sn11)/float(n11i))*((_sn_1-_sn11)/float(n11i+_sn-_sn1_-_sn_1)) )
                _sn11 = n11i
                return _log_sprob
            if n11i==_sn11-1:
                _log_sprob = _log_sprob + \
                        log ( ((_sn11)/float(_sn1_-n11i))*((_sn11+_sn-_sn1_-_sn_1)/float(_sn_1-n11i)) )
                _sn11 = n11i
                return _log_sprob
        _sn11 = n11i
    else:
        _sn11 = n11i
        _sn1_=n1_i
        _sn_1=n_1i
        _sn=ni
    _log_sprob = log_hyper_323(_sn11,_sn1_,_sn_1,_sn)
    return _log_sprob

def log_hyper(n11):
    return log_hyper0(n11,0,0,0)

def tgtmocnt(infile, version, tgtmodict):
    ifile = open(infile, 'rt')
    reader = csv.reader(ifile, delimiter = '\t')
    ifile.seek(0)
    for row in reader:
	if version in row[0]:
	    ## create a dictionary with (target id, motif id) as keys and motif counts as values 
	    tgtmodict[(row[1],row[0].split('_')[0])] += 1
	#else:
	    #print 'wrong version of motif database :/'
    return tgtmodict

def mocnt(tgtmodict, modict):
    while tgtmodict:
	## create a dictionary of hits per motif for a set of sequences
	(tgt, mo), cnt = tgtmodict.popitem()
	if cnt > 0:
	    modict[mo] += 1
    return modict	
    
def main(argv):
    if len(argv) < 4:
        sys.stderr.write("Usage: %s exp_fimo_output_file bkgd_fimo_output_file motif_list_file\n" % argv[0])
        return 1
    if not os.path.isfile(argv[1]):
        sys.stderr.write('Error: experiment_fimo_output_file %r was not found!\n' % argv[1])
        return 1
    if not os.path.isfile(argv[2]):
        sys.stderr.write('Error: background_fimo_output_file %r was not found!\n' % argv[2])
        return 1
    if not os.path.isfile(argv[3]):
        sys.stderr.write('Error: motif_list_file %r was not found!\n' % argv[3])
        return 1
    
    infile_exp = sys.argv[1]
    infile_bkgd = sys.argv[2]
    infile_motiflist = sys.argv[3]

    (path,fname) = os.path.split(infile_exp)
    (n,ext) = os.path.splitext(fname)

    ifile_motiflist = open(infile_motiflist,'rt')
    reader_motiflist = csv.reader(ifile_motiflist, delimiter = '\t')
       
    ofile = open(os.path.join(path, 'test_out'), 'wt')
    writer = csv.writer(ofile, delimiter = '\t')

    ##(path1,fname1) = os.path.split(infile1)
    ##(path2,fname2) = os.path.split(infile2)
    
    mlist = []
    for row in reader_motiflist:
        motifid = row[0].split(' ')[-1].split('_')[0]
        if not motifid in mlist:
            mlist.append(motifid)
 
    modictbkgd = defaultdict(int)
    modictexp = defaultdict(int)
    tgtmodictbkgd = defaultdict(int)
    tgtmodictexp = defaultdict(int)
    tgtmocountbkgd = tgtmocnt(infile_bkgd, '_0.80', tgtmodictbkgd)
    #with open('tgtmocountbkgd.txt','w') as f:
	#pickle.dump(tgtmocountbkgd,f)
    #with open('tgtmocountbkgd.txt','r') as f:
	#tgtmocountbkgd = pickle.load(f)
    tgtmocountexp = tgtmocnt(infile_exp, '_0.80', tgtmodictexp)
    #with open('tgtmocountexp.txt','w') as f:
        #pickle.dump(tgtmocountexp,f)
    #with open('tgtmocountexp.txt','r') as f:
        #tgtmocountexp = pickle.load(f)
    #print tgtmocountbkgd

    modictbkgd = mocnt(tgtmocountbkgd, modictbkgd)   
    #with open('modictbkgd.txt','w') as f:
        #pickle.dump(modictbkgd,f)
    with open('motif2tfname.txt','r') as f:
        mydict1 = pickle.load(f)
    print len(mydict1)
    modictexp = mocnt(tgtmocountexp, modictexp)
    #for mo in modictbkgd.keys():
        #print len(modictbkgd[mo])
    #expmocount = defaultdict(int)
    #ifile_exp.seek(0)
    #for row in reader_exp:
	#if '_0.62' in row[0]:
	#if row[1] == 'Tcrb':
	    #expmocount[row[0]] += 1
    #print expmocount
    

    pval = defaultdict(float)
    log_pvalue = defaultdict(float)
       
    P = 3730 ##total number of exp sequences
    N = 68654#59249#60100 ##total number of bkgd sequences

    for mo in mlist:
	#print mo 
	if mo in modictexp.keys():
	    p = modictexp[mo]	
	else:
	    p = 0
	if mo in modictbkgd.keys():
	    n = modictbkgd[mo]
	else:
	    n = 0
	## calculate empirical p-val
	#pval[mo] = (counter[mo]+1.0)/(len(modictbkgd[mo])+1)    
	## calculate hypergeometric p-val
	#log_pvalue[mo] = getLogFETPvalue(p, P, n, N, -5)
        log_pvalue[mo] = log_getFETprob(N-n, n, P-p, p)[4];
	# print the pvalue
	if isinstance(log_pvalue[mo],tuple):
            pvalue = sprint_logx(log_pvalue[mo][1], 3, "%6.3fe%-5.0f")
	else:
	    pvalue = sprint_logx(log_pvalue[mo], 3, "%6.3fe%-5.0f")
        print >> sys.stdout, log_pvalue[mo], pvalue, p, P, n, N, mo, mydict1[mo]
	newrow = [mo, mydict1[mo], pvalue, p, P, n, N]
	writer.writerows([newrow])
    	
    ofile.close()
	
if __name__=='__main__':
    sys.exit(main(sys.argv))




