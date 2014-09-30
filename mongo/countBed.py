import os, sys, csv
from array import array
from copy import copy
from collections import defaultdict
import time
#from ngs_plumbing.intervals import Interval, IntervalList
import struct

def compressBed4(bedFile, ctName, bbFile):
	"""compress bedFiles into binary"""
        with open(bbFile, 'wb') as f:
                bedFile.seek(0)
                dataString = bedFile.read()
                st = struct.Struct('cciif')
		dataList = dataString.split('\n')
                for d in dataList[:-1]:
			chrom = d.split('\t')[0]
                        if len(chrom) == 4:
                                chrom = '0'+chrom[-1]
                        elif len(chrom) == 5:
                                chrom = chrom[3:5]
			#"0-based" bedGraph
			start, end = int(d.split('\t')[1]), int(d.split('\t')[2])
			score = float(d.split('\t')[3])
                        ls = list(chrom)
                        line = ls + [start,end,score]
                        l = st.pack(*line)
                        f.write(l)
        return 0

def getBinBedCoord(bbFile, expName):
	"""get chromosomal coordinates and values stored in dictionaries
		unpack from binary file"""
	coordDict = defaultdict(lambda: defaultdict(list))
	valuesDict = defaultdict(lambda: defaultdict(list))
	st = struct.Struct('cciif')
	stSize = struct.calcsize('cciif')
	with open(bbFile, 'rb') as f:
		flag = 1
		while flag:
			chunk = f.read(stSize)
			if stSize > len(chunk):
				flag = 0
			else:
				data = st.unpack(chunk)
				if data[0] == '0':
					chrom = 'chr'+data[1]
				else:
					chrom = 'chr'+data[0]+data[1]
				start,end,value = data[2]+1, data[3]+1, data[4]
				coordDict[expName][chrom].extend(xrange(start,end))
				valuesDict[expName][chrom].extend([value]*(end-start))
	return coordDict, valuesDict

def getBed4Coord(bedFile, expName):
        """get chromosomal coordinates and values 
		in BED4 format stored in dictionaries"""
        coordDict = defaultdict(lambda : defaultdict(list))
        valuesDict = defaultdict(lambda : defaultdict(list))
        for row in bedFile:
                #print row
		if len(row) == 4:##assert row format
                	chrom = row[0]
			start, end = int(row[1])+1, int(row[2])+1
			value = float(row[3])
			coordDict[expName][chrom].extend(xrange(start,end))
			valuesDict[expName][chrom].extend([value]*(end-start))
		else:
			print "unmatched file format"
        return coordDict, valuesDict

def getBed6Coord(bedFile, expName):
	"""get chromosomal coordinates and values
		in BED6	format stored in dictionaries"""
	coordDict = defaultdict(lambda : defaultdict(list))
	valuesDict = defaultdict(lambda: defaultdict(list))
	for row in bedFile:
		if len(row) == 6:##assert format
			chrom = row[0]
			start, end = int(row[1])+1, int(row[2])+1
			reg = row[3]
			value = float(row[4])
			#info = row[-1]
			coordDict[expName][chrom].extend(xrange(start,end))
			valuesDict[expName][chrom].extend([value]*(end-start))
		else:
			print "unmatched file format"
	return coordDict, valuesDict

def getBed6Anno(bedFile, expName):
	"""get chromosomal interval information
		in BED6 format stored in dictionaries"""
	annoIntvlDict = defaultdict(lambda : defaultdict(tuple))
	for row in bedFile:
		if len(row) == 6:
			chrom = row[0]
			anno = row[3]
			start, end = int(row[1])+1, int(row[2])
			#start, end = int(row[1]),int(row[2])
			value = float(row[4])
			annoIntvlDict[chrom][(start,end)] = (anno, value)
			#annoIntvlDict[chrom][Interval(start,end)] = (anno, value)
		else:
			print "unmatched file format"
	return annoIntvlDict

def sortInterval(annoIntvlDict):
    """sort all intervals on a chromosome"""
    ##deprecated-- requires the use of interval object
    intervalDict = {}
    for chrom in annoIntvlDict:
        itl = IntervalList( intv for intv in annoIntvlDict[chrom])
        itl.sort()
        intervalDict[chrom] = itl
    return intervalDict

def sortStart(annoIntvlDict):
    intervalDict = defaultdict(list)
    for chrom in annoIntvlDict:
	itl = list((start,end) for (start,end) in annoIntvlDict[chrom])
	#print itl
	intervalDict[chrom] = sorted(itl,key=lambda x: x[0])
    return intervalDict

def sortEnd(annoIntvlDict):
    intervalDict = defaultdict(list)
    for chrom in annoIntvlDict:
	itl = list((start,end) for (start,end) in annoIntvlDict[chrom])
	intervalDict[chrom] = sorted(itl,key=lambda x: x[1])
    return intervalDict

def getMotifAnno(annoIntvlDict,intervalStartDict,intervalEndDict,motifChrom,motifStart,motifEnd,window):
    """push excluded regions from BED6 into a list
	if motif falls into it"""
    regionList = []
    valueList = []
    #itl = IntervalList( intv for intv in geneRangeDict[motifChrom])
    try:
	startList = [start for (start,end) in intervalStartDict[motifChrom]]
	endList = [end for (start,end) in intervalEndDict[motifChrom]]
	intervalStartList = [(start,end) for (start,end) in intervalStartDict[motifChrom]]
	intervalEndList = [(start,end) for (start,end) in intervalEndDict[motifChrom]]
	iStart = bisect(startList, motifEnd)
	iEnd = bisect(endList,motifStart)
	s = intervalStartList[:iStart]
	e = intervalEndList[iEnd:]
	overlapping = list(set(s)&set(e))
        #intervals = intervalDict[motifChrom]##for short list of Bed6 Anno, brute force is Ok
        #motifInterval = Interval(motifStart, motifEnd)
        #overlapping = [ x for x in intervals \
	#	if ((x[1]+window)>motifInterval[0] and (x[0]-window)<motifInterval[0]) \
	#	or ((x[1]+window)>motifInterval[1] and (x[0]-window)<motifInterval[1]) ]

        #print 'overlapping',overlapping
        if not len(overlapping) == 0:
            for x in overlapping:
                (anno, value) = annoIntvlDict[motifChrom][x]
                if not anno in regionList:
                    regionList.append(anno)
                if not (anno, value) in valueList:
                    valueList.append((anno,value))
        #if len(targetList) > 1:
            #print targetList, motifChrom, motifInterval
    except KeyError:
        pass
    return regionList, valueList

def buildBedHist(chrom,coordDict,valuesDict,ctName):
        """build histogram from bed input with discontinuous interval values
		same as building histograms from variableStep wig files
                e.g.
          		chr1	3000000	3000059	1
			chr1	3000059	3000060	0.1
			chr1	3000060	3000061	0.0625
			chr1	3000061	3000062	0.0384615
			chr1	3000062	3000068	0.05
			chr1	3000068	3000069	0.0149254
			chr1	3000069	3000071	0.0188679
			chr1	3000071	3000073	0.030303
			chr1	3000073	3000077	0.0238095
			chr1	3000077	3000079	0.030303
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
	n = len(xs) # last ending elemented is appended
	ll, rr = 0, n-1
	while ll <= rr:
		m = (ll+rr)/2
		if xs[m] > start:
			rr = m - 1
		else:
			ll = m + 1
	li = rr#ll	# get the left side index
	ll, rr = 0, n-1
	while ll <= rr:
		m = (ll+rr)/2
		if xs[m] > end:
			rr = m - 1
		else:
			ll = m + 1
	ri = rr # get the right side index
	if ri < li:
		return 0.0	# nothing in between, maybe wrong [start, end]?
	sum = sums[ri-1] - sums[li-1]
	sum -= (start - xs[li]) * xvals[li]		# remove the left extra area
	sum += (end - xs[ri] + 1) * xvals[ri]	# remove the right extra area

	##brute force
	check = False
	if check:
		ans = 0
		for i in range(0,n):
			if start <= xs[i] and (i==n-1 or end >= xs[i+1]):
				if i!=n-1:
					ans += (xs[i+1]-xs[i]) * xvals[i]
			elif xs[i] <= start and (i==n-1 or start < xs[i+1]):
				if i!=n-1:
					ans += (xs[i+1] - start) * xvals[i]
			elif xs[i] <= end and (i==n-1 or end < xs[i+1]):
				ans += (end-xs[i]+1) * xvals[i]
		if ans != sum:
			print "No - Wrong Answer:", ans, sum
			return -1
		else:
			print "Yes", ans

	size = end - start + 1
	avg = sum/size

	return avg, size
	
def main(argv):
	if len(argv) < 3:
		sys.stderr.write("Usage: %s gff_file wig_file \n" % argv[0])
		return 1
	if not os.path.isfile(argv[1]):
		sys.stderr.write('Error: gff_file %r was not found!\n' % argv[1])
		return 1
	if not os.path.isfile(argv[2]):
		sys.stderr.write('Error: wig_file %r was not found!\n' % argv[2])
		return 1
	
	# example usage
	#time1 = time.time()
	wig = open(sys.argv[2],'rt')
	#wigFile = csv.reader(wig, delimiter = '\t')
	ctName = "mapability"
	bwFile = '/home/xc406/data/mongodbtest/testbedGraph.bb'
	#coordDict, valuesDict = getBedCoord(wigFile, ctName)
	#dataType = 'phyloP46wayPrimate'
	#compressBed4(wig, ctName, bwFile)
	getBinBed4Coord(bwFile, ctName)
	#stepDict, startDict, valuesDict = getFixStart(wigFile, dataType)#'phyloP46wayPrimate')
	#dataType = 'phyloP46wayPrimate'
	chrom = 'chr1'
	#print chrom
	#arrayDict = buildFixHist(chrom, stepDict, startDict, valuesDict, dataType)
	#start = startDict[dataType][chrom]
	#for i in xrange(len(start)):
	#	xs, xvals, sums = arrayDict[start[i]]
	#	print "xs", xs
	#	print "xvals", xvals
	#	print "sums", sums
	#	avg, size = queryHist(xs, xvals, sums, 60084, 60090)
		#avg, size = queryHist(xs, xvals, sums, 61338037, 61338046)
	#	print avg, size
	#	break
	#print 'wig processing time',time.time()-time1

	#time2 = time.time()
#	gff = open(sys.argv[1],'rt')
	##gffFile = HTSeq.GFF_Reader(gff)
#	gffFile = csv.reader(gff, delimiter = '\t')
#	features = getRange(gffFile)
#	intvlen = len(features)
	#print 'gff processing time',time.time()-time2
	
	#time3 = time.time()	
	##build the arrays (3 needed)
	#arrayDict=defaultdict(list)
	#arrayDict[chrom] = buildHist(chrom,coordDict,valuesDict,ctName)
	#xs, xvals, sums = arrayDict[chrom]
	#start, end = 3000055, 3000070 ##include both boudaries
	#avg, size = queryHist(xs, xvals, sums, start, end)
	#print avg, size

#	for i in xrange(intvlen):
#		chrom,start,end = features[i][0],features[i][1],features[i][2]
#		if not chrom in arrayDict:
#			arrayDict[chrom] = buildHist(chrom,coordDict,valuesDict)
#		xs, xvals, sums = arrayDict[chrom]
		#print chrom,'index', xs, len(xs)
		#print 'values', xvals, len(xvals)
		#print 'sums',sums,len(sums)
#		avg = queryHist(xs, xvals, sums, start, end)[0]					# pass the three arrays to the query
		#avg2 = queryHist(xs, xvals, sums, -100, 108)			# query again	
		#print 'average count:',chrom,start,end,avg
		#print 'count time', time.time()-time3
#		print avg
	return 0

if __name__=='__main__':
    sys.exit(main(sys.argv))
