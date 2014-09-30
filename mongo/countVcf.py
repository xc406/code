import os, sys, csv, re, gzip
from array import array
from copy import copy
from collections import defaultdict
import time
import struct
import itertools
from operator import mul
from fractions import Fraction

def nCk(n,k):
	return int(reduce(mul,(Fraction(n-i,i+1) for i in range(k)),1))

def pack(f,chrnum,numLine,num):
	header = struct.Struct('cci')
	numst = struct.Struct('if'*numLine)
	if len(chrnum)==1:
	    chrom = '0'+chrnum
	elif len(chrnum)==2:
	    chrom = chrnum
	ls = list(chrom)
	ls.append(numLine)
	h = header.pack(*ls)
	f.write(h)
	nl = numst.pack(*num)
	f.write(nl)
	return 0


def compressVcf(vcfFile, expName, bwFile, bwFileSnp, bwFileIndel):
	"""encode vcf files in binary format as header(cci)+numString(if)"""
	"""where cci=chrnum+numLine ii=pos+p##+ac+an+alt"""
	#header = struc.Struct('cci')
	#numst = struct.Struct('iii')
	with open(bwFile,'wb') as f, open(bwFileSnp, 'wb') as fSnp, open(bwFileIndel, 'wb') as fIndel:
	    d = gzip.open(vcfFile)
	    pchrnum = '0'
	    while True:
		l = d.readline()
		if l=='':##check for EOF
		    #write the remaining data out
		    print "pack last chromosome", pchrnum
		    pack(f,pchrnum,numLine,num)
		    pack(fSnp,pchrnum,numLineSnp,numSnp)
		    pack(fIndel,pchrnum,numLineIndel,numIndel)
		    break
		elif not l.startswith('#'):
		    chrnum = l.split('\t')[0]
		    if not pchrnum == chrnum:
			if not pchrnum == '0':
			    print "pack", pchrnum
			    pack(f,pchrnum,numLine,num)
			    pack(fSnp,pchrnum,numLineSnp,numSnp)
			    pack(fIndel,pchrnum,numLineIndel,numIndel)
			##write out data for one chrom and reset	
			numLine, numLineSnp, numLineIndel = 0,0,0
			num, numSnp, numIndel  = [], [], []
			pchrnum = chrnum
		    ra = l.split('\t')[3]
		    aa = l.split('\t')[4].split(',')
		    v = 'Snp'
		    if len(ra) > 1:
			v = 'Indel'
		    else:
		    	for a in aa:
			    if len(a) > 1:
			    	v = 'Indel'
		    pos = int(l.split('\t')[1])
		    an = int(l.split('AN=')[-1].split(';')[0])
		    numLine += 1
		    try:
		    	ac = int(l.split('AC=')[-1].split(';')[0])
			#ac = min(ac,an-ac)##ensure minor allele freq
			p = nCk(ac, 1)*nCk(an-ac,1)/float(nCk(an,2))
		    except ValueError:#if multiple snps
			aclist = l.split('AC=')[-1].split(';')[0].split(',')
			#for i in xrange(n):##calculate the upper triangular mat
			acs = [int(x) for x in aclist]
			acs.append(an-sum(acs))
			n = len(acs)
			for ac1,ac2 in itertools.combinations(acs,2):
			    p += nCk(ac1,1)*nCk(ac2,1)
			p /=float(nCk(an,n))
		    num.extend([pos,p])
		    if v == "Snp":
			numSnp.extend([pos,p])
			numLineSnp += 1
		    elif v == "Indel":
			numIndel.extend([pos,p])
			numLineIndel += 1
	return 0 

def getBinVarCoord(bwFile, expName):
        """get chromosomal coordinates and values stored in dictionaries
                unpack from bw files """
        coordDict = defaultdict(lambda: defaultdict(list))
        valuesDict = defaultdict(lambda: defaultdict(list))
        numst = struct.Struct('if')
        chunkSize = struct.calcsize('if')
        headerst = struct.Struct('cci')
        headerSize = struct.calcsize('cci')
        with open(bwFile,'rb') as f:
                header = f.read(headerSize)
                h1,h2,s = headerst.unpack(header)
                if h1 == '0':
                        chrom = 'chr'+h2
                else:
                        chrom = 'chr'+h1+h2
                while True:
                        for i in xrange(s):
                                chunk = f.read(chunkSize)
                                pos, p = numst.unpack(chunk)
				coordDict[expName][chrom].append(pos)
                                valuesDict[expName][chrom].append(p)
                        header = f.read(headerSize)
                        l = len(header)
                        if l < headerSize:
                                break
                        else:
                                h1, h2, s = headerst.unpack(header)
                                if h1 == '0':
                                        chrom = 'chr'+h2
                                else:
                                        chrom = 'chr'+h1+h2
        return coordDict, valuesDict

def buildVarHist(chrom,coordDict,valuesDict,ctName):#coord,values):
        """build histogram from wig input with discontinuous values
                e.g.
                variableStep chrom=chr13
                19021446        1.00
                19022345        1.00
                19022949        3.00
                19022956        1.00
                19025166        5.00
                19025399        1.00
                19025986        1.00
                19026391        1.00
                19026727        1.00
                """
        coord = [-1] + coordDict[ctName][chrom] # prepend left boundary
        values = [0] + valuesDict[ctName][chrom]
        n = len(values)
        x, lastx, lastval = coord[0], coord[0]-1, -100000000
        xs, xvals = array("i"), array("d")      # clear and reset histogram, index is integer, value is double
        for i in xrange(n):             # get the histogram
                x = coord[i]
                val = values[i]
                if x != lastx+1:
                        xs.append(lastx+1)
                        xvals.append(0)
                        #lastx = lastx+1
                        lastval = 0
                if val != lastval:
                        xs.append(x)
                        xvals.append(val)
                lastx = x
                lastval = val
        xs.append(coord[n-1]+1) # append right boundary
        xvals.append(0)
        nlen = len(xs)
        sums = [0] * nlen
        sums[0] = xvals[0]*(xs[1]-xs[0])
        for i in xrange(1,nlen-1):
                sums[i] = sums[i-1] + (xs[i+1]-xs[i])*xvals[i]
        sums[nlen-1] = sums[nlen-2]
        #print "xs",xs
        #print "xvals",xvals
        #print "sums",sums
        return xs, xvals, sums

def queryHist(xs, xvals, sums, start, end):
        """query histogram to get average value of a defined genomic region"""
        #if start < xs[0]:
        #       print 'start out of range'
        #elif end > xs[-1]:
        #       print 'end out of range'
        n = len(xs) # last ending elemented is appended
        ll, rr = 0, n-1
        while ll <= rr:
                m = (ll+rr)/2
                if xs[m] > start:
                        rr = m - 1
                else:
                        ll = m + 1
        li = rr#ll      # get the left side index
        ll, rr = 0, n-1
        while ll <= rr:
                m = (ll+rr)/2
                if xs[m] > end:
                        rr = m - 1
                else:
                        ll = m + 1
        ri = rr # get the right side index
        if ri < li:
                return 0.0      # nothing in between, maybe wrong [start, end]?
        sum = sums[ri-1] - sums[li-1]
        sum -= (start - xs[li]) * xvals[li]             # remove the left extra area
        sum += (end - xs[ri] + 1) * xvals[ri]   # remove the right extra area

        size = end - start + 1
	avg = sum/size

        return avg, size, sum

def main(argv):
	if len(argv) < 3:
		sys.stderr.write("Usage: %s vcf_file bw_file \n" % argv[0])
		return 1
	if not os.path.isfile(argv[1]):
		sys.stderr.write('Error: vcf_file %r was not found!\n' % argv[1])
		return 1
	if not os.path.exists(argv[2]):
		sys.stderr.write('Error: path bw_file %r was not found!\n' % argv[2])
		return 1

	vcfFile = sys.argv[1]
	(vcfPath,vcfFname) = os.path.split(sys.argv[2])
	bwFile = os.path.join(sys.argv[2],'test.vcf.bw')
	expName = "1000genomes"
	startTime = time.clock()
	compressVcf(vcfFile, expName, bwFile)
	print time.clock()-startTime
	coordDict, valueDict = getBinVarCoord(bwFile, expName)		
	print time.clock()-startTime
	chrom, start,end= 'chr1', 233050,233095
	if not chrom in arrayDict:
		arrayDict[chrom] = buildVarHist(chrom, coordDict,valuesDict, expName)
	xs, xvals, sums = arrayDict[chrom]
	s1,s2 = queryHist(xs,xvals,sums,start,end)
	s3,s4 = queryHist(xs,xvals,sums,536865,536896)

	print s1,s2,'chr1', 233050,233095
	print s3,s4, 536865,536896
	print time.clock()-startTime

if __name__=='__main__':
    sys.exit(main(sys.argv))
