from pymongo import MongoClient
from pymongo.errors import DuplicateKeyError
import os, csv, sys, glob, gzip, time
from pymongo import ASCENDING, DESCENDING
import numpy as np
from bisect import bisect
from collections import defaultdict
import countWig, countBed, countVcf, writeMongo
from Bio import SeqIO
from scipy.stats import binom

def updateCount(path, tfName):
    """update count data in discontinuous variableStep wiggle format"""
    for infile in glob.glob(os.path.join(path,"*.wig")):
        #(wigpath,wigfilename) = os.path.split(infile)
	(wigfilename,ext) = os.path.splitext(infile)
        ##depends on the data type and source
        expName = "Dnase"##or add an expName parser line
        ctName = wigfilename.split('EncodeUwDnase')[-1].split('Aln')[0]
        wigFile = open(infile,'rt')
        #wig = csv.reader(wigFile,delimiter='\t')
	countWig.compressVarWig(wigFile, expName, wigfilename)
        coordDict, valuesDict = countWig.getBinVarCoord(wigfilename+motifChrom+'.bw',ctName)
        arrayDict = defaultdict(list)
        cursor = mcollection.find({"tf_name": tfName})
        for test in cursor:
            motifChrom = test["motif_genomic_regions_info"]["chr"]
	    motifStart = test["motif_genomic_regions_info"]["start"] 
	    motifEnd = test["motif_genomic_regions_info"]["end"]
            if not motifChrom in arrayDict:
                arrayDict[motifChrom] = countWig.buildHist(motifChrom,coordDict,valuesDict,ctName)
            xs, xvals, sums = arrayDict[motifChrom]
            count = countWig.queryHist(xs, xvals, sums, motifStart, motifEnd)[0]
            #print count
	    #mcollection.update({"_id":test["_id"]},{"$set":{expName: count}}, upsert = True)
            test["ct_info"]["accessibility_score"][expName] = count
            mcollection.save(test)
    return 0

def updateFOS(path, tfName, motifChrom, method="NSD", flankWin=35):
    """calculate fos from discontinuous variableStep wiggle files
	with two method options:
		NSD/Binomial test"""
    for infile in glob.glob(os.path.join(path,"*.wig")):
	#(wigpath,wigfile) = os.path.split(infile)
	(wigfilename,ext) = os.path.splitext(infile)
	expName = "FOS"
	ctName = wigfilename.split('EncodeUwDgf')[-1].split('Aln')[0]
	#wigFile = open(infile,'rt')
	#wig = csv.reader(wigFile,delimiter='\t')
	#bwFile = os.path.join(path,wigfilename+'.bw')
	#countWig.compressVarWig(wigFile, expName, wigfilename)
	coordDict, valuesDict = countWig.getBinVarCoord(wigfilename+motifChrom+'.bw',ctName)
	arrayDict = defaultdict(list)
	cursor = mcollection.find({"tf_name": tfName,
		"motif_score":{"$lt":1e-4},
		"motif_genomic_regions_info.chr": motifChrom})
	for test in cursor:
	    if not motifChrom in arrayDict:
		arrayDict[motifChrom] = countWig.buildVarHist(motifChrom,coordDict,valuesDict,ctName)
	    xs, xvals, sums = arrayDict[motifChrom]
	    motifStart = test["motif_genomic_regions_info"]["start"]
	    motifEnd = test["motif_genomic_regions_info"]["end"]
	    flankL = max(0, motifStart - flankWin)
	    flankR = motifEnd + flankWin
	    countTotL = countWig.queryHist(xs, xvals, sums, flankL, motifEnd)[2]
	    countTotR = countWig.queryHist(xs, xvals, sums, motifStart, flankR)[2]
	    countCent = countWig.queryHist(xs, xvals, sums, motifStart, motifEnd)[2]
	    count = countWig.queryHist(xs, xvals, sums, motifStart-100, motifEnd+100)[2]
	    if method == "NSD":
		try:
		    fos = np.sqrt((countTot-countCent)/countTot)-np.sqrt(countCent/countTot)
		except ZeroDivisionError:
		    fos = 0 
	    elif method == "Binom":
		try:
		    fos = min(1 - binom.cdf(countCent,countTotL,float(motifEnd-motifStart)/(motifEnd-flankL)), 1 - binom.cdf(countCent,countTotR,float(motifEnd-motifStart)/(flankR-motifStart)))
		except ZeroDivisionError:
		    fos = 0
	    #mcollection.update({"_id":test["_id"]},{"$set":{"motif_ct_info.fos": {ctName:{method:fos}}}}, upsert = True)
	    if fos > 0.95 and count-countCent > 18:#(flankR-flankL-(motifEnd-motifStart)):
	    	print motifChrom+'\t'+str(motifStart)+'\t'+str(motifEnd)+'\t'+str(fos)
    return 0 	

def updatePi(path, tfName, motifChrom):
    """calculate nucleotide diversity for a given genomic region"""
    cursor = mcollection.find({"tf_name": tfName, "motif_genomic_regions_info.chr": motifChrom})
    for infile in glob.glob(os.path.join(path,"*v5.20130502.sites.vcf.gz")):
    	(vcfpath,vcffile) = os.path.split(infile)
	(vcffilename,vcfext) = os.path.splitext(vcffile)

	bwFile = os.path.join(path,vcffilename+'.bw')
	bwFileSnp = os.path.join(path,vcffilename+'_snp.bw')
	bwFileIndel = os.path.join(path,vcffilename+'_indel.bw')

	expName = '1000genomes'

	#countVcf.compressVcf(infile, expName, bwFile, bwFileSnp, bwFileIndel)

	#coordDict, valuesDict = countVcf.getBinVarCoord(bwFile, expName)
	#arrayDict = defaultdict(list)

	#coordDictSnp, valuesDictSnp = countVcf.getBinVarCoord(bwFileSnp, expName)
	#arrayDictSnp = defaultdict(list)

	coordDictIndel, valuesDictIndel = countVcf.getBinVarCoord(bwFileIndel, expName)
	arrayDictIndel = defaultdict(list)

    	for test in cursor:
	    #motifChrom = test["motif_genomic_regions_info"]["chr"]
	    motifStart = test["motif_genomic_regions_info"]["start"]
	    motifEnd = test["motif_genomic_regions_info"]["end"]
	    #if not motifChrom in arrayDict:
                #arrayDict[motifChrom] = countVcf.buildVarHist(motifChrom, coordDict, valuesDict, expName)
	    #if not motifChrom in arrayDictSnp:
		#arrayDictSnp[motifChrom] = countVcf.buildVarHist(motifChrom, coordDictSnp, valuesDictSnp, expName)
	    if not motifChrom in arrayDictIndel:
		arrayDictIndel[motifChrom] = countVcf.buildVarHist(motifChrom, coordDictIndel, valuesDictIndel, expName)
	    #xs, xvals, sums = arrayDict[motifChrom]
	    #avg = countVcf.queryHist(xs,xvals,sums,motifStart,motifEnd)[0]
	    #xs, xvals, sums = arrayDictSnp[motifChrom]
	    #avgSnp = countVcf.queryHist(xs,xvals,sums,motifStart,motifEnd)[0]
	    xs, xvals, sums = arrayDictIndel[motifChrom]
	    avgIndel = countVcf.queryHist(xs,xvals,sums,motifStart,motifEnd)[0]
	    mcollection.update({"_id":test["_id"]},{"$set":{#"motif_mapability_info.piTot": avg, 
			#"motif_mapability_info.piSnp": avgSnp, 
			"motif_mapability_info.piIndel": avgIndel}}, upsert = True)
	    print "pi=", avgIndel#,motifEnd-motifStart
    return 0

def updateGC(path, tfName, motifChrom):
    """update GC content for a given genomic region"""
    cursor = mcollection.find({"tf_name": tfName,"motif_genomic_regions_info.chr": motifChrom,
		"motif_score":{"$lt":1e-3}})
    try:
	chrnum = int(motifChrom[3:])
	fname = 'mm9_chrom%02d' % chrnum
    except:
	fname = 'mm9_chrom20'
    with open(os.path.join(path,fname)) as f:
	seq_record = SeqIO.read(f,"fasta")
        for test in cursor:
	    #motifChrom = test["motif_genomic_regions_info"]["chr"]
	    motifStart = test["motif_genomic_regions_info"]["start"]
	    motifEnd = test["motif_genomic_regions_info"]["end"]
	#with open(os.path.join(path,motifChrom+'.fa.masked')) as f:
	    #seq_record = SeqIO.read(f,"fasta")
	    seq = seq_record.seq[motifStart-1:motifEnd]##(]
	    size = float(motifEnd-motifStart+1)

	    if not size == float(seq.count('C')+seq.count('G')+seq.count('A')+seq.count('T')):
		gc = (seq.count('C')+seq.count('G')+seq.count('c')+seq.count('g'))/size
		if size == seq.count('N'):
		    print 'N repeat found', seq
		    mcollection.update({"_id":test["_id"]},{"$set":{"motif_mapability_info.gc_content": ('NA','repeat')}}, upsert = True)
		if size == float(seq.count('c')+seq.count('g')+seq.count('a')+seq.count('t')):
	            #print 'repeat found', gc, seq
		    mcollection.update({"_id":test["_id"]},{"$set":{"motif_mapability_info.gc_content": (gc,'repeat')}}, upsert = True)
		if size > seq.count('N') and seq.count('N') > 0:
		    print 'partial repeat found', seq
		    mcollection.update({"_id":test["_id"]},{"$set":{"motif_mapability_info.gc_content": ('NA','partial repeat')}}, upsert = True)
		elif seq.count('N') == 0:
		    #print 'partial repeat found', gc, seq
		    mcollection.update({"_id":test["_id"]},{"$set":{"motif_mapability_info.gc_content": (gc,'partial repeat')}}, upsert = True)
	    else:
		gc = (seq.count('C')+seq.count('G'))/size
		#print 'not repeats', gc, seq
		mcollection.update({"_id":test["_id"]},{"$set":{"motif_mapability_info.gc_content": (gc,'non-repeat')}}, upsert = True)
    return 0

def updateCons(path, tfName):
    """update conservation scores in gzipped fixedStep wiggle format"""
    for infile in glob.glob(os.path.join(path, "*.gz")):
	(wigpath,wigfilename) = os.path.split(infile)
	chrom = wigfilename.split('.')[0]
	consName = '_'.join(wigfilename.split('.')[1:3])
	#print chrom, tfName, consName
	with gzip.open(infile) as wigFile:
	    wig = csv.reader(wigFile,delimiter='\t')
	    stepDict, startDict, valuesDict = countWig.getFixStart(wig,consName)
	    start = startDict[consName][chrom]
	    arrayDict = countWig.buildFixHist(chrom,stepDict,startDict,valuesDict,consName)
	    cursor = mcollection.find({"tf_name": tfName, 
				"motif_genomic_regions_info.chr": chrom})
	    for test in cursor:
	        motifStart, motifEnd = test["motif_genomic_regions_info"]["start"], test["motif_genomic_regions_info"]["end"]
	    	avg = 0
		#print motifStart, motifEnd
		startlist = [start[i] for i in xrange(len(start)-1) if (motifStart >= start[i] and motifStart < start[i+1]) or (motifEnd >= start[i] and motifEnd < start[i+1])]
		if motifEnd > start[-1]:
		    startlist.append(start[-1])
	    	for i in xrange(len(startlist)):
		 #   if avg != 0 and motifEnd < start[i]:
		  #  	break ##fall into range and break out
		    if avg != 0:
			if motifEnd >= startlist[i]:##cases of partial overlap need to renormalize over two fragments
		    	    xs, xvals, sums = arrayDict[startlist[i]]
		    	    avg = avg * (startlist[i] - motifStart) + (countWig.queryHist(xs, 
				xvals, sums, motifStart, motifEnd)[0] * (motifEnd - startlist[i] + 1)) / (motifEnd - motifStart + 1)	
			else:
			    break
		    elif avg == 0:
		   	xs, xvals, sums = arrayDict[startlist[i]]
		   	avg = countWig.queryHist(xs, xvals, sums, motifStart, motifEnd)[0]
	    	#if avg > 0:
	            #print avg, motifStart, motifEnd
	    	mcollection.update({"_id":test["_id"]},{"$set":{"motif_cons_info":{consName: avg}}}, upsert = True)
		#else:
		    #mcollection.update({"_id":test["_id"]},{"$set":{"motif_cons_info":{consName: avg}}}, upsert = True)
    return 0

def updateMap(path, tfName, window):
    """update mapability scores from bedGraph files"""
    for infile in glob.glob(os.path.join(path, "*.bedGraph")):##better in bedGraph
	(mappath,mapfilename) = os.path.split(infile)
	expName = mapfilename.split('.')[0]
	#print expName
	with open(infile,'rt') as bedGraphFile:
	    bedGraph = csv.reader(bedGraphFile, delimiter = '\t')
	    coordDict, valuesDict = countBed.getBed4Coord(bedGraph, expName)
	    arrayDict=defaultdict(list)
	    cursor = mcollection.find({"tf_name":tfName})
	    for test in cursor:
		motifChrom, motifStart, motifEnd = test["motif_genomic_regions_info"]["chr"], \
                                                test["motif_genomic_regions_info"]["start"], \
                                                test["motif_genomic_regions_info"]["end"]
		if not motifChrom in arrayDict:
		    arrayDict[motifChrom] = countBed.buildBedHist(motifChrom, coordDict, valuesDict, expName)
		xs, xvals, sums = arrayDict[motifChrom]
		avg, size = countBed.queryHist(xs, xvals, sums, motifStart-window, motifEnd+window)
		
		if avg != 0:
		    print motifChrom, motifStart, motifEnd
		    print avg
 		mcollection.update({"_id":test["_id"]},{"$set":{"motif_mapability_info":{"score":{expName:avg}}}}, upsert = True)
    return 0

def updateExcludedRegions(path,tfName,window):
    """mark motifs if it overlaps ENCODE excluded regions"""
    for infile in glob.glob(os.path.join(path,"*.bed.gz")):
	(regpath,regfilename) = os.path.split(infile)
	expName = regfilename.split('.')[0]
	with gzip.open(infile,'rt') as bedFile:
	    bed = csv.reader(bedFile, delimiter = '\t') 
	    annoIntvlDict = countBed.getBed6Anno(bed,expName)
	    intervalDict = countBed.sortInterval(annoIntvlDict)
	    cursor = mcollection.find({"tf_name": tfName})
	    for test in cursor:
		motifChrom, motifStart, motifEnd = test["motif_genomic_regions_info"]["chr"], \
						test["motif_genomic_regions_info"]["start"], \
						test["motif_genomic_regions_info"]["end"]
	        regionList, valueList = countBed.getMotifAnno(annoIntvlDict,
				intervalDict,motifChrom,motifStart,motifEnd,window)
		if regionList != []:
			print regionList, valueList, motifChrom, motifStart, motifEnd
		mcollection.update({"_id": test["_id"]}, 
			{"$set": {"motif_mapability_info":{"exclude": regionList}}}, upsert = True)
	    	
    return 0

def updateMotifGenomicRegions(path):
    """clean overlapping motif entries and update motif genomic loci"""
    for infile in glob.glob(os.path.join(path,"mm9_wg_1e3_90_*14.txt")):
        (fimopath,fimofilename) = os.path.split(infile)
        tfName = fimofilename.split('mm9_wg_1e3_90_')[-1].split('fimoout')[0].split('_')[0]
        f = open(infile,'rt')
        fimo = csv.reader(f,delimiter='\t')
        pchrom = 'chr0'
        rowlist = []
        for row in fimo:
            if not "pattern" in row[0]:
                ##fimo output is 1-based, so don't have to change anything in mongodb or gff
                chrom, start, end, strand, pval = row[1], int(row[2]), int(row[3]), row[4], float(row[6])
                size = end - start + 1
                if chrom == pchrom:
                ##########strand specific############ 
		#if control for strand-specificity, split fimo files as + and -, and run separately
                    if pend in xrange(start, end+1) or start in xrange(pstart, pend+1):
                        if abs(pend-start) > min(size, psize)/2.0:
                           rowlist.append(row)
                        else:
                            if len(rowlist) > 1:
                                wrow = rowlist[minP(rowlist)]
                                ##best of the overlapping motifs
				##current line slightly overlapping but ok
                            elif len(rowlist) == 1:
                                wrow = prow
                                ##previous line non-overlapping
				##current line slightly overlapping but ok
                            rowlist = []
                            rowlist.append(row)
                            wmid, wchrom, wstart, wend, wstrand, wpval = wrow[0], wrow[1], int(wrow[2]), int(wrow[3]), wrow[4], float(wrow[6])
                            m = mcollection.find_one({"motif_id": wmid, "tf_name": tfName})
                            mcollection.insert({"motif_id": wmid, 
						"tf_name": tfName, 
						#"motif_gene_mapping_info": m["motif_gene_mapping_info"],
						#"motif_ct_info": m["motif_ct_info"], 
						#"motif_tf_info": m["motif_tf_info"],
						#"motif_mapability_info": m["motif_mapability_info"],
						#"motif_cons_info": m["motif_cons_info"], 
						"motif_genomic_regions_info": 
						{"chr": wchrom, 
						"start": wstart, 
						"end": wend, 
						"strand": wstrand}, 
						"motif_score": wpval})
                            mcollection.remove({"tf_name": tfName, "motif_score": None})##remove original placeholder entry
                            #ofile.writelines('\t'.join(modRow(wrow))+'\n')
                    else:
                        if len(rowlist) > 1:
                            wrow = rowlist[minP(rowlist)]
                            ##best of the overlapping motifs
			    ##current line non-overlapping
                        elif len(rowlist) == 1:
                            wrow = prow
                            ##previous and current line non-overlapping
                        rowlist = []##reset container
                        rowlist.append(row)
                        wmid, wchrom, wstart, wend, wstrand, wpval = wrow[0], wrow[1], int(wrow[2]), int(wrow[3]), wrow[4], float(wrow[6])
                        m = mcollection.find_one({"motif_id": wmid, "tf_name": tfName})
                        mcollection.insert({"motif_id": wmid, 
					    "tf_name": tfName, 
					    #"motif_gene_mapping_info": m["motif_gene_mapping_info"],
					    #"motif_ct_info": m["motif_ct_info"], 
					    #"motif_tf_info": m["motif_tf_info"], 
					    #"motif_mapability_info": m["motif_mapability_info"],
					    #"motif_cons_info": m["motif_cons_info"],
					    "motif_genomic_regions_info": 
					     {"chr": wchrom, 
					      "start": wstart, 
					      "end": wend, 
					      "strand": wstrand}, 
					    "motif_score": wpval})
                prow = row
                psize = size
                if len(rowlist) == 0:
                    rowlist.append(prow)##push first row in, during looping this should never be empty
                pchrom, pstart, pend, pstrand, ppval = chrom, start, end, strand, pval

def getRefSeqDict(refSeqReader):
    """store transcripts loci info in dictionaries"""
    tssDict = defaultdict(list)
    geneNameDict = defaultdict(tuple)
    sizeDict = defaultdict(int)
    geneRangeDict = defaultdict(lambda : defaultdict(tuple))
    for row in refSeqReader:
	try:
	    chrom, strand, txStart, txEnd, geneName, transcriptId = row[2], row[3], int(row[4]), int(row[5]), row[12], row[1]
	    size = txEnd - txStart
	    if strand == '+':
	        if not geneName in sizeDict:
	    	    tssDict[chrom].append(txStart)
	    	    geneNameDict[(chrom,txStart)] = (geneName, transcriptId)
		    sizeDict[geneName] = size
		    #geneRangeDict[chrom][Interval(txStart,txEnd)] = (geneName, transcriptId)
		    geneRangeDict[chrom][(txStart,txEnd)] = (geneName, transcriptId)
	        else:
		    if size > sizeDict[geneName]:
		    	tssDict[chrom].append(txStart)
		    	geneNameDict[(chrom,txStart)] = (geneName, transcriptId)
		    	sizeDict[geneName] = size
			geneRangeDict[chrom][(txStart,txEnd)] = (geneName, transcriptId)
	    else:
	    	if not geneName in sizeDict:
		    tssDict[chrom].append(txEnd)
		    geneNameDict[(chrom,txEnd)] = (geneName, transcriptId)
		    sizeDict[geneName] = size
		    geneRangeDict[chrom][(txStart,txEnd)] = (geneName, transcriptId)
	    	else:
		    if size > sizeDict[geneName]:
	    	    	tssDict[chrom].append(txEnd)
	    	    	geneNameDict[(chrom,txEnd)] = (geneName, transcriptId)
		    	sizeDict[geneName] = size
			geneRangeDict[chrom][(txStart,txEnd)] = (geneName, transcriptId)
	except ValueError:
	    pass
    return tssDict, geneNameDict, geneRangeDict

def getClosestGene(tssDict,geneNameDict,motifChrom,motifStart,motifEnd):
    """find closest gene feature to motif entries"""
    for chrom in tssDict:
    	tssDict[chrom].sort()##don't sort beforehand, change the original dictionary
    motifMid = motifStart+(motifEnd-motifStart)/2
    try:
    	i = bisect(tssDict[motifChrom],motifMid)
    	if i < len(tssDict[motifChrom]) and i > 0:
            if abs(motifMid - tssDict[motifChrom][i-1]) >= abs(tssDict[motifChrom][i] - motifMid):
            	(geneName, transcriptId) = geneNameDict[(motifChrom,tssDict[motifChrom][i])]
            	dist = abs(tssDict[motifChrom][i] - motifMid)
            else:
            	(geneName, transcriptId) = geneNameDict[(motifChrom,tssDict[motifChrom][i-1])]
            	dist = abs(motifMid - tssDict[motifChrom][i-1])
    	elif i == len(tssDict[motifChrom]):
            if not i == 0:
            	(geneName, transcriptId) = geneNameDict[(motifChrom,tssDict[motifChrom][i-1])]
            	dist = abs(motifMid - tssDict[motifChrom][i-1])
            else:
            	(geneName, transcriptId) = (None, None)
            	dist = None
    	elif i == 0:
            (geneName, transcriptId) = geneNameDict[(motifChrom,tssDict[motifChrom][i])]
            dist = abs(motifMid - tssDict[motifChrom][i])
    except KeyError:
	geneName, dist, transcriptId = None, None, None
    return geneName, dist, transcriptId

def getTargetGene(geneRangeDict, intervalStartDict, intervalEndDict, motifChrom, motifStart, motifEnd, txWindow):
    """push motif target genes into a list"""
    targetList = []
    transcriptList = []
    #itl = IntervalList( intv for intv in geneRangeDict[motifChrom])
    #itl.sort()
    try:
    	startList = [start-1-txWindow for (start,end) in intervalStartDict[motifChrom]]
	##-1 changes 1-based encoding into 0-based
	endList = [end+txWindow for (start,end) in intervalEndDict[motifChrom]]
	intervalStartList = [(start,end) for (start,end) in intervalStartDict[motifChrom]]
	intervalEndList = [(start,end) for (start,end) in intervalEndDict[motifChrom]]
	iStart = bisect(startList, motifEnd)
	iEnd = bisect(endList, motifStart)
	s = intervalStartList[:iStart]##bed intervals are 0-based, so overlapping with bed start is not accounted for
	e = intervalEndList[iEnd:]##motifStart is 1-based, so overlapping with bed end is acccounted for
	overlapping = list(set(s)&set(e))
    	#motifInterval = Interval(motifStart, motifEnd)
    	#overlapping = [ x for x in intervals if (x[1]>motifInterval[0] and x[0]<motifInterval[0]) \
	#or (x[1]>motifInterval[1] and x[0]<motifInterval[1]) ]
        #print 'overlapping',overlapping
	
    	if not len(overlapping) == 0:
	    for x in overlapping:
		#print x
	    	(geneName, transcriptId) = geneRangeDict[motifChrom][x]
	    	if not geneName in targetList:
	            targetList.append(geneName)
		if not (geneName,transcriptId) in transcriptList:
		    transcriptList.append((geneName,transcriptId))
	#if len(targetList) > 1:
	    #print targetList, motifChrom, motifInterval
    except KeyError:
	pass
    return targetList, transcriptList

def modRow(row):
    """fix motif size in gff file-- not sure what it does now"""
    mlen = str(int(row[-1].split(' ')[1].split('_')[3])+2)
    group = [row[0],row[3],row[4],mlen,row[6]]
    row[-1] = 'gene_id ' + '_'.join(group) #+ '; sequence ' + prow[-1].split(' ')[-1]
    return row

def minP(rowlist):
    """return the index of the row with the lowest pval"""
    pval = float(rowlist[0][6])
    minp = 0
    for p in xrange(len(rowlist)):
	if pval > float(rowlist[p][6]):
	    minp = p
	    pval = float(rowlist[p][6])
    return minp

def main(argv):
    if len(argv) < 3:
        sys.stderr.write("Usage: %s motif_tf_info_file path-to-fimo-output\n" % argv[0])
        return 1
    if not os.path.isfile(argv[1]):
        sys.stderr.write('Error: motif_info_file %r was not found!\n' % argv[1])
        return 1
    if not os.path.exists(argv[2]):
        sys.stderr.write('Error: path-to-fimo-output %r was not found!\n' % argv[2])
        return 1 

    server = 'localhost'
    port = 27017
    client = MongoClient(server, port)
    db = client["mm9"]
    global mcollection 
    mcollection = db["motif_instance_hughes_test"]
    ##drop collection
    #c["mm9"].drop_collection("motif_instance_hughes_test")
    #mcollection.remove()
    #db.drop_collection('motif_instance_hughes_test')
    #mcollection = db["motif_instance_hughes_test"]
    #print 'clean', mcollection.count()
    #index collections
#    mcollection.ensure_index("motif_id",name="m_id",unique=False,background=True)
#    mcollection.ensure_index("tf_name",name="tf_name",unique=False,background=True)
    #collection.ensure_index("motif_type",name="motif_type",unique=True,drop_dups=True,background=True)
#    mcollection.ensure_index("genomic_regions_gene_mapping.genelist10kb",name = "target_gene",unique=False,background=True)
#    mcollection.ensure_index("genomic_regions_gene_mapping.closest_gene",name = "closest_gene",unique=False,background=True)
#    mcollection.ensure_index("motif_score", name = "motif_score", unique = False, background = True)
#    mcollection.ensure_index("motif_tf_info.motif_type", name = "motif_type", unique = False, background = True)
#    mcollection.ensure_index("motif_tf_info.msource_type", name = "msource_type", unique = False, background = True)
#    mcollection.ensure_index("motif_tf_info.tf_status", name = "tf_status", unique = False, background = True)
#    mcollection.ensure_index("motif_tf_info.msource_id", name = "project_name", unique = False, background = True)
#    mcollection.ensure_index("motif_tf_info.family_name", name = "family_name", unique = False, background = True)
#    mcollection.ensure_index([("tf_name", ASCENDING),
#			      ("motif_gene_mapping_info.genelist10kb", DESCENDING)],
#				name="network_edge", unique=False, background=True)
    #index genomic regions
    #mcollection.ensure_index([("motif_genomic_regions_info.chr",DESCENDING),
#				("motif_genomic_regions_info.start",DESCENDING),
				#("motif_genomic_regions_info.end",ASCENDING),
#				("motif_genomic_regions_info.strand",DESCENDING)],
#				name="genomic_regions",unique=False,background=True)
    #print collection
    #print 'done indexing'
    infile = sys.argv[1]#'/home/xc406/data/hg19motifs90/TF_Information90hg19.txt'
    path = sys.argv[2]

    ifile = open(infile,'rt')
    tf_info = csv.reader(ifile, delimiter = '\t')
    #mlist = []
    for row in tf_info:
	try:
		dbd_count = int(row[11])
	except ValueError:
		dbd_count = None
	try:
		msource_year = int(row[18])
	except ValueError:
		msource_year = None
	#print row
	motif_instance = {
		"motif_id": row[3],
		"tf_name": row[6],
		"motif_score": None,
		"motif_tf_info":{	
			"species_name": row[7],
			"tf_status": row[8], #direct or indirect
			"family_name": row[9],
			"dbds": row[10],
			"dbd_count": dbd_count,
			"dbid": row[12],
			"motif_type": row[14],
			"msource_id": row[15],
			"msource_type": row[16],
			"msource_author": row[17],
			"msource_year": msource_year,
			"pmid": row[19] #citation
		},
		"motif_genomic_regions_info":{
			"chr": None,
			"start": None,
			"end": None,
			"strand": None
		},
#		"motif_mapability_info":{
#				"exclude": [],
#				"score": [],
#				"gc_content": None
#		},
#		"motif_cons_info":{
#				"phylop_euarchontoglires": None,
#				"phylop_mammals": None,
#				"phylop_vertebrate": None,
#				"phastCons_euarchontoglires": None,
#				"phastCons_mammals": None,
#				"phastCons_vertebrate": None
				#"SNP_diversity": None,
				#"Indel_diversity": None
#		},
#		"motif_gene_mapping_info":{
#			"closest_gene": None,
#			"feature": None,#intergenic or 3' 5'
#			"dist_tss": None,
#			"genelist10kb": [],##list of gene--center of motif fall in gene plus and minus 10kb
#			"transcriptidlist10kb": []
			#"epu_id": #boolean		
#		},
#		"motif_ct_info":{
#			"ct_name": None,
#			"ct_type": None, #normal/cancerous/cellline/primarycell
#			"accessibility_score": {}, #{type: log likelihood score}
#			#"accessibility_type": [], #dhs, dgf, faire
#			"chip_score": {}, #tf: pval(overlapping peaks)
#			"h3k4me3_score": None,
#			"h3k4me1_score": None,
#			"h3k27ac_score": None,
#			"p300_score": None,
#			"pol2_score": None
#		},	
	}

	#mlist.append(motif_instance)
    	#print len(mlist)

 #       try:	
		#print motif_instance	
#		mcollection.insert(motif_instance)
#		del motif_instance
		#print 'inserted' #collection
 #   	except DuplicateKeyError:
		#print 'dup'
#		pass

    #cursor = mcollection.find()
    #print 'before entering genomic region info ', mcollection.count()
    #c = iter(cursor)

    ##clean overlapping motif entries
    startTime = time.clock()
#    updateMotifGenomicRegions(path)
    
    #getCount(wigpath,"Hes5")

    ##write gff
#    cursor = mcollection.find({"tf_name": "Stat3",
#			"motif_genomic_regions_info.chr":{"$nin": ["chrM", "chr13_random","chr16_random","chr17_random",
#			"chr1_random", "chr3_random", "chr4_random", "chr5_random", "chr7_random",
#			"chr8_random", "chr9_random", "chrUn_random", "chrX_random", "chrY_random"]}})
#    print 'updated count ', cursor.count(), 'total count after update', mcollection.count(), 'update time', time.clock() - startTime

    #updateGC(path,"Stat3")
    #i, o = getPi('/data/1000genomes/ALL.autosomes.phase3_shapeit2_mvncall_integrated_v4.20130502.sites.vcf.gz')
    #updatePi('/data/1000genomes/','Stat3','chr6')
    updateFOS(path, 'Ctcf', 'chr17', method="Binom")
#    ofile = open('/home/xc406/data/mongodbtest/Runx1.gff','wt')
#    gffWriter = csv.writer(ofile, delimiter='\t')

    ##update gene features
#    refSeqFile = open('/home/xc406/data/mm9_refseq_June_2014.txt','rt')
#    refSeqReader = csv.reader(refSeqFile, delimiter='\t')
#    tssDict, geneNameDict, geneRangeDict = getRefSeqDict(refSeqReader)
    #cursor = mcollection.find({"tf_name":"Stat3"})
    #intervalDict = sortInterval(geneRangeDict)

#    intervalStartDict = countBed.sortStart(geneRangeDict)
#    intervalEndDict = countBed.sortEnd(geneRangeDict)
    #with client.start_request():##open update
    #print "Zic1 motif: {0}".format(test)
#    for test in cursor:
#        motifChrom, motifStart, motifEnd = test["motif_genomic_regions_info"]["chr"], test["motif_genomic_regions_info"]["start"], test["motif_genomic_regions_info"]["end"]
    	#print closestGene(tssDict,geneNameDict,motifChrom,motifStart,motifEnd)[0][0]
    	#print closestGene(tssDict,geneNameDict,motifChrom,motifStart,motifEnd)[1]
    	#print test["motif_id"], test["motif_genomic_regions_info"]["chr"],test["motif_genomic_regions_info"]["start"]
	#startTime1 = time.clock()
#	t = getTargetGene(geneRangeDict,intervalStartDict,intervalEndDict,motifChrom, motifStart, motifEnd, 10000)
	#endTime1 = time.clock()
#	#print "all target mapping", t, motifChrom, motifStart, motifEnd, endTime1 - startTime1
#	closest = getClosestGene(tssDict,geneNameDict,motifChrom,motifStart,motifEnd)
	#endTime2 = time.time()
	#print "closest gene mapping time", closest, endTime2 - endTime1 
    	#mcollection.update(test,{"$set":{"genomic_regions_gene_mapping":{"closest_gene": closestGene(tssDict,geneNameDict,motifChrom,motifStart,motifEnd)[0][0], 
#			"dist_tss": closestGene(tssDict,geneNameDict,motifChrom,motifStart,motifEnd)[1]}}})
	#if len(closestGene(tssDict,geneNameDict,motifChrom,motifStart,motifEnd)[0]) > 1:
	    #print closestGene(tssDict,geneNameDict,motifChrom,motifStart,motifEnd)[0], motifChrom, motifStart, motifEnd
#	mcollection.update({"_id":test["_id"]},{"$set":{"genomic_regions_gene_mapping":{"closest_gene": (closest[0],closest[2]), 
#			"dist_tss": closest[1], "genelist10kb": t[0], "transcriptidlist10kb": t[1]}}}, upsert = True)
    	#test["genomic_regions_gene_mapping"]["closest_gene"] = (closest[0],closest[2])#closestGene(tssDict,geneNameDict,motifChrom,motifStart,motifEnd)[0][0] 
   	#test["genomic_regions_gene_mapping"]["dist_tss"] = closest[1]#closestGene(tssDict,geneNameDict,motifChrom,motifStart,motifEnd)[1]
	#test["genomic_regions_gene_mapping"]["genelist10kb"] = t[0]
	#test["genomic_regions_gene_mapping"]["transcriptidlist10kb"] = t[1]#getTargetGene(geneRangeDict,geneNameDict,motifChrom, motifStart, motifEnd, 0)
    	#mcollection.save(test)
	#print "Hes5 motif: {0}".format(test)	

    #updateCons(path,"Hes5")
    #updateExcludedRegions(path,"Hes5",0)
    #updateMap(path, "Hes5", 0)
    #print 'update time', time.time() - startTime2
    #cursor = mcollection.find({"tf_name":"Runx1",
#			"motif_genomic_regions_info.chr":{"$nin": ["chrM", "chr13_random","chr16_random","chr17_random",
#			"chr1_random", "chr3_random", "chr4_random", "chr5_random", "chr7_random",
#			"chr8_random", "chr9_random", "chrUn_random", "chrX_random", "chrY_random"]}})
#    writeMongo.makeGff(cursor,gffWriter,0)
#    ofile.close()
    #testupdate = mcollection.find_one({"tf_name":"Zic1"})#{"motif_id": test["motif_id"], 
		#"motif_genomic_regions_info":{"chr": test["motif_genomic_regions_info"]["chr"], "start": test["motif_genomic_regions_info"]["start"]}})
    #print "Zic1 motif: {0}".format(test)
    #print "Zic1 motif update: {0}".format(testupdate)
    print 'total time', time.clock() - startTime
    #mcollection.insert({motif_id: item['motif_id']})
    #c = cursor.next()
    #print cursor.next()
    ##update motif_genomic_regions_info
    #db.motif_instance_hughes_test.insert({"motif_tf_info": db.motif_instance_hughes_test.findOne({'tf_name':'Irf4'}).motif_tf_info})

if __name__=='__main__':
    sys.exit(main(sys.argv))
