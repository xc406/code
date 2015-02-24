from pymongo import MongoClient
from pymongo.errors import DuplicateKeyError
import re, os, csv, sys, glob, gzip, time, bz2
from pymongo import ASCENDING, DESCENDING
import numpy as np
from bisect import bisect
from collections import defaultdict
import countWig, countBed, countVcf, writeMongo
from Bio import SeqIO
from scipy.stats import binom

def updateCount(path, db, motifChrom='chr17', window=100):
    """update count data in discontinuous variableStep wiggle format"""
    mcollection = db["hg19"+motifChrom]
    for infile in glob.glob(os.path.join(path,"wgEncodeUwDgf*_cut"+motifChrom)):
        #(wigpath,wigfilename) = os.path.split(infile)
	#(wigfilename,ext) = os.path.splitext(infile)
	wigfilename = infile.split(motifChrom)[0]
        ##depends on the data type and source
        expName = "dgf"##or add an expName parser line
        ctName = wigfilename.split('EncodeUwDgf')[-1].split('_cut')[0]
        wigFile = open(infile,'rt')
        #wig = csv.reader(wigFile,delimiter='\t')
	if not os.path.isfile(wigfilename+motifChrom+'.bw'):
	    countWig.compressVarWig(wigFile, expName, wigfilename)
        coordDict, valuesDict = countWig.getBinVarCoord(wigfilename+motifChrom+'.bw',ctName)
        arrayDict = defaultdict(list)
        cursor = mcollection.find({"tf_name":{"$in": ["IRF3","MAFK","NFYA","SIN3A","ZNF384"]}})
        for test in cursor:
            #motifChrom = test["genomic_region"]["chr"]
	    motifStart = test["genomic_region"]["start"] 
	    motifEnd = test["genomic_region"]["end"]
            if not motifChrom in arrayDict:
                arrayDict[motifChrom] = countWig.buildVarHist(motifChrom,coordDict,valuesDict,ctName)
            xs, xvals, sums = arrayDict[motifChrom]
	    #print motifStart, motifEnd, len(xs)
            count = countWig.queryHist(xs, xvals, sums, motifStart-window, motifEnd+window, varWindow=True)[0]
            #print count
	    mcollection.update({"_id":test["_id"]},{"$set":{expName+"."+ctName: count}}, upsert = True)
            #test["ct_info"]["accessibility_score"][expName] = count
            #mcollection.save(test)
    return 0

def updateFOS(path, db, motifChrom="chr17", dgfCutoff=11, method="Binom", stranded=False):#, flankWin=35):
    """calculate fos from discontinuous variableStep wiggle files
	with two method options:
		NSD/Binomial test"""
    mcollection = db["hg19"+motifChrom]
    print 'updating fos', motifChrom
    expName="fos"
    for infile in glob.glob(os.path.join(path,"wgEncodeUwDnaseGm12878Aln_cut"+motifChrom)):
	#(wigpath,wigfile) = os.path.split(infile)
	#(wigfilename,ext) = os.path.splitext(infile)
	wigfilename = re.split(motifChrom,infile)[0]#"_._cut",infile)[0]
	ctName = wigfilename.split('EncodeUwDnase')[-1].split('Aln')[0]
	if not stranded:
		wigFile = open(infile,'rt')
		#wig = csv.reader(wigFile,delimiter='\t')
		#bwFile = os.path.join(path,wigfilename+'.bw')
		#countWig.compressVarWig(wigFile, expName, wigfilename)
		bwFile = wigfilename+motifChrom+'.bw'
		if not os.path.isfile(bwFile):
	    		countWig.compressVarWig(wigFile, expName, wigfilename)
		coordDict, valuesDict = countWig.getBinVarCoord(bwFile,ctName)
		arrayDict = defaultdict(list)
		cursor = mcollection.find()#{"tf_name":{"$in": ["IRF3","MAFK","NFYA","SIN3A","ZNF384"]}})
		for test in cursor:
	    		if not motifChrom in arrayDict:
				arrayDict[motifChrom] = countWig.buildVarHist(motifChrom,coordDict,valuesDict,ctName)
	    		xs, xvals, sums = arrayDict[motifChrom]
	    		motifStart = test["genomic_region"]["start"]
	    		motifEnd = test["genomic_region"]["end"]
	    		flankWin = round((motifEnd - motifStart + 1)*1.75)
	    		flankL = max(0, motifStart - flankWin)
	    		flankR = motifEnd + flankWin
	    		countTotL = countWig.queryHist(xs, xvals, sums, flankL, motifEnd, varWindow=True)[2]
	    		countTotR = countWig.queryHist(xs, xvals, sums, motifStart, flankR, varWindow=True)[2]
	    		countCent = countWig.queryHist(xs, xvals, sums, motifStart, motifEnd)[2]
	    		count = countWig.queryHist(xs, xvals, sums, motifStart-100, motifEnd+100, varWindow=True)[2]
	    		if method == "NSD":
				try:
		    			fos = np.sqrt((count-countCent)/count)-np.sqrt(countCent/count)
				except ZeroDivisionError:
		    			fos = 0 
	    		elif method == "Binom":
				try:
		    			fos = min(1 - binom.cdf(countCent,countTotL,float(motifEnd-motifStart)/(motifEnd-flankL)), 
					1 - binom.cdf(countCent,countTotR,float(motifEnd-motifStart)/(flankR-motifStart)))
				except ZeroDivisionError:
		    			fos = 0
	    		if fos > 0.95 and count-countCent > dgfCutoff:#(flankR-flankL-(motifEnd-motifStart)):
				mcollection.update({"_id":test["_id"]},{"$set":{"dnase.fos": fos}}, upsert = True)
	    			#print motifChrom+'\t'+str(motifStart)+'\t'+str(motifEnd)+'\t'+str(fos)
    	else:
		infilep = os.path.join(path,"wgEncodeUwDnaseGm12878Aln_p_cut.wig")
		infilen = os.path.join(path,"wgEncodeUwDnaseGm12878Aln_n_cut.wig")
		wigFilep = open(infilep,'rt')
		wigFilen = open(infilen,'rt')
		bwFilep = wigfilename+'_p_cut'+motifChrom+'.bw'
		bwFilen = wigfilename+'_n_cut'+motifChrom+'.bw'
		if not os.path.isfile(bwFilep):
			countWig.compressVarWig(wigFilep,expName,wigfilename+'_p_cut')
		if not os.path.isfile(bwFilen):
			countWig.compressVarWig(wigFilen,expName,wigfilename+'_n_cut')
		coordDictp, valueDictp = countWig.getBinVarCoord(bwFilep,ctName)
		coordDictn, valueDictn = countWig.getBinVarCoord(bwFilen,ctName)
		arrayDictp = defaultdict(list)
		arrayDictn = defaultdict(list)
		cursor = mcollection.find()
		for test in cursor:
			if not motifChrom in arrayDictp:
				arrayDictp[motifChrom] = countWig.buildVarHist(motifChrom,coordDictp,valueDictp,ctName)
			if not motifChrom in arrayDictn:
				arrayDictn[motifChrom] = countWig.buildVarHist(motifChrom,coordDictn,valueDictn,ctName)
			xs, xvals, sums = arrayDictp[motifChrom]
			motifStart,motifEnd = test["genomic_region"]["start"],test["genomic_region"]["end"]
			flankWin = round((motifEnd - motifStart + 1)*1.75)
			flankL= max(0, motifStart-flankWin)
			flankR = motifEnd + flankWin
			countTotLp = countWig.queryHist(xs, xvals, sums, flankL, motifEnd, varWindow=True)[2]
			countCentp = countWig.queryHist(xs, xvals, sums, motifStart, motifEnd)[2]
			countp = countWig.queryHist(xs, xvals, sums, motifStart-100, motifEnd+100, varWindow=True)[2]
			xs, xvals, sums = arrayDictn[motifChrom] 
			countTotRn = countWig.queryHist(xs, xvals, sums, motifStart, flankR, varWindow=True)[2]
			countCentn = countWig.queryHist(xs, xvals, sums, motifStart, motifEnd)[2]
			countn = countWig.queryHist(xs, xvals, sums, motifStart-100, motifEnd+100, varWindow=True)[2]
			if method == "Binom":
				try:
					pp = 1 - binom.cdf(countCentp,countTotLp,float(motifEnd-motifStart)/(motifEnd-flankL))
					pn = 1 - binom.cdf(countCentn,countTotRn,float(motifEnd-motifStart)/(flankR-motifStart))
					fos = pp*pn
				except ZeroDivisionError:
					fos = 0
			if fos > 0.95 and countp+countn-countCentp-countCentn > dgfCutoff:
				mcollection.update({"_id":test["_id"]},{"$set":{"dnase.fosStrand": fos}}, upsert = True)
				print motifChrom+'\t'+str(motifStart)+'\t'+str(motifEnd)+'\t'+str(fos)
    return 0 	

def updatePi(path, db, motifChrom="chr17"):
    """calculate nucleotide diversity for a given genomic region"""
    mcollection = db["hg19"+motifChrom]
    print 'updating pi'
    #mcollection.ensure_index("pi",name="nucleotide_diversity",unique=False,background=True)
    cursor = mcollection.find({"tf_name":{"$in": ["IRF3","MAFK","NFYA","SIN3A","ZNF384"]}})
	#{"tf_name": tfName, "motif_genomic_regions_info.chr": motifChrom})
    for infile in glob.glob(os.path.join(path,"*v5.20130502.sites.vcf.bw")):
    	(vcfpath,vcffile) = os.path.split(infile)
	(vcffilename,vcfext) = os.path.splitext(vcffile)

	#bwFile = os.path.join(path,vcffilename+'.bw')
	bwFileSnp = os.path.join(path,vcffilename+'_snp.bw')
	bwFileIndel = os.path.join(path,vcffilename+'_indel.bw')

	expName = '1000genomes'

	#countVcf.compressVcf(infile, expName, bwFile, bwFileSnp, bwFileIndel)

#	coordDict, valuesDict = countVcf.getBinVarCoord(bwFile, expName)
#	arrayDict = defaultdict(list)

	coordDictSnp, valuesDictSnp = countVcf.getBinVarCoord(bwFileSnp, expName)
	arrayDictSnp = defaultdict(list)

	coordDictIndel, valuesDictIndel = countVcf.getBinVarCoord(bwFileIndel, expName)
	arrayDictIndel = defaultdict(list)

    	for test in cursor:
	    #motifChrom = test["motif_genomic_regions_info"]["chr"]
	    motifStart = test["genomic_region"]["start"]
	    motifEnd = test["genomic_region"]["end"]
#	    if not motifChrom in arrayDict:
 #               arrayDict[motifChrom] = countVcf.buildVarHist(motifChrom, coordDict, valuesDict, expName)
	    if not motifChrom in arrayDictSnp:
		arrayDictSnp[motifChrom] = countVcf.buildVarHist(motifChrom, coordDictSnp, valuesDictSnp, expName)
	    if not motifChrom in arrayDictIndel:
		arrayDictIndel[motifChrom] = countVcf.buildVarHist(motifChrom, coordDictIndel, valuesDictIndel, expName)
#	    xs, xvals, sums = arrayDict[motifChrom]
#	    avg = countVcf.queryHist(xs,xvals,sums,motifStart,motifEnd)[0]
	    xs, xvals, sums = arrayDictSnp[motifChrom]
	    avgSnp = countVcf.queryHist(xs,xvals,sums,motifStart,motifEnd)[0]
	    xs, xvals, sums = arrayDictIndel[motifChrom]
	    avgIndel = countVcf.queryHist(xs,xvals,sums,motifStart,motifEnd)[0]
	    mcollection.update({"_id":test["_id"]},{"$set":{#"map.pi": avg, 
			"cons.piSnp": avgSnp, 
			"cons.piIndel": avgIndel}}, upsert = True)
	    #print "pi=", avgIndel#,motifEnd-motifStart
    return 0

def updateGC(path, db, motifChrom='chr17'):
    """update GC content for a given genomic region
	(path-to-fasta-files db collection-name)"""
    mcollection = db["hg19"+motifChrom]
    mcollection.ensure_index("gc",name="gc_content",unique=False,background=True)
    print 'updating GC'
    cursor = mcollection.find({"tf_name":{"$in": ["IRF3","MAFK","NFYA","SIN3A","ZNF384"]}})
	#{"tf_name": tfName,"motif_genomic_regions_info.chr": motifChrom,
		#"motif_score":{"$lt":1e-3}})
    try:
	chrnum = int(motifChrom[3:])
	fname = 'hg19_chrom%02d' % chrnum
    except:
	fname = 'hg19_chrom23'
    with open(os.path.join(path,fname)) as f:
	seq_record = SeqIO.read(f,"fasta")
        for test in cursor:
	    #motifChrom = test["motif_genomic_regions_info"]["chr"]
	    motifStart = test["genomic_region"]["start"]
	    motifEnd = test["genomic_region"]["end"]
	    seq = seq_record.seq[motifStart-1:motifEnd]##(]
	    size = float(motifEnd-motifStart+1)

	    if not size == float(seq.count('C')+seq.count('G')+seq.count('A')+seq.count('T')):
		gc = (seq.count('C')+seq.count('G')+seq.count('c')+seq.count('g'))/size
		if size == seq.count('N'):
		    print 'N repeat found', seq
		    #mcollection.update({"_id":test["_id"]},{"$set":{"motif_mapability_info.gc_content": ('NA','repeat')}}, upsert = True)
		if size == float(seq.count('c')+seq.count('g')+seq.count('a')+seq.count('t')):
	            #print 'repeat found', gc, seq
		    mcollection.update({"_id":test["_id"]},{"$set":{"gc": (gc,'r')}}, upsert = True)
		if size > seq.count('N') and seq.count('N') > 0:
		    print 'partial repeat found', seq
		    #mcollection.update({"_id":test["_id"]},{"$set":{"motif_mapability_info.gc_content": ('NA','partial repeat')}}, upsert = True)
		elif seq.count('N') == 0:
		    #print 'partial repeat found', gc, seq
		    mcollection.update({"_id":test["_id"]},{"$set":{"gc": (gc,'pr')}}, upsert = True)
	    else:
		gc = (seq.count('C')+seq.count('G'))/size
		#print 'not repeats', gc, seq
		mcollection.update({"_id":test["_id"]},{"$set":{"gc": (gc,'nr')}}, upsert = True)
    return 0

def updateCons(path, db, motifChrom='chr17'):
    """update conservation scores in gzipped fixedStep wiggle format"""
    mcollection = db["hg19"+motifChrom]
    for infile in glob.glob(os.path.join(path, motifChrom+"*.wigFix.gz")):
	(wigpath,wigfilename) = os.path.split(infile)
	chrom = wigfilename.split('.')[0]
	consName = '_'.join(wigfilename.split('.')[1:-2])
	print 'updating', consName
	#print chrom, tfName, consName
	with gzip.open(infile) as wigFile:
	    #wig = csv.reader(wigFile,delimiter='\t')
	    bwFile = os.path.join(wigpath,motifChrom+"."+consName+'.bw')
	    if not os.path.isfile(bwFile):
		countWig.compressFixWig(wigFile, consName, bwFile)
	    stepDict, startDict, valuesDict = countWig.getBinFixStart(bwFile,consName)
	    start = startDict[consName][chrom]
	    arrayDict = countWig.buildFixHist(chrom,stepDict,startDict,valuesDict,consName)
	    cursor = mcollection.find({"tf_name":{"$in": ["IRF3","MAFK","NFYA","SIN3A","ZNF384"]}})
			#{"tf_name": tfName, 
				#"genomic_region.chr": chrom})
	    print mcollection.count()
	    #num = 0
	    #avg = 0
	    for test in cursor:
		motifStart, motifEnd = test["genomic_region"]["start"], test["genomic_region"]["end"]
		#num += 1
		#print avg#motifStart, motifEnd, num	
		avg = 0
		startlist = [start[i] for i in xrange(len(start)-1) if (motifStart >= start[i] and motifStart < start[i+1]) or (motifEnd >= start[i] and motifEnd < start[i+1])]
		if motifEnd > start[-1]:
		    startlist.append(start[-1])
		for i in xrange(len(startlist)):
		    #if avg != 0:
			#if motifEnd >= startlist[i]:##cases of partial overlap need to renormalize over two fragments
		    ss = startlist[i]
		    xs, xvals, sums, ll = arrayDict[ss]
		    if motifStart < ss <= motifEnd <= ss+ll-1:##left out, right in
			if avg == 'NA' and i == len(startlist)-1:
			    avg = 0
			avg += countWig.queryHist(xs,xvals, sums, ss, motifEnd)[0] *(motifEnd - ss + 1) /(motifEnd - motifStart + 1)
		    elif ss <= motifStart < motifEnd <= ss+ll-1:##in array
			avg = countWig.queryHist(xs,xvals, sums, motifStart, motifEnd)[0]
		    elif motifStart < ss and ss+ll-1 < motifEnd:##motif > array
			if avg == 'NA':
			    avg = 0 
			avg += countWig.queryHist(xs,xvals, sums, ss, ss+ll-1)[0] * ll /(motifEnd - motifStart + 1)
		    elif ss <= motifStart <= ss+ll-1 < motifEnd:##left in, right out
			if avg == 'NA' and i == len(startlist)-1:
			    avg = 0
			avg += countWig.queryHist(xs,xvals, sums, motifStart, ss+ll-1)[0] *(ss + ll - motifStart) /(motifEnd - motifStart + 1)
		    elif ss+ll-1 < motifStart:
			if avg == 0:
			    #print 'here-->', motifStart, motifEnd, ss, ll
			    avg = 'NA'
		    elif motifEnd < ss:
			print "this should not happen-- motifStart < motifEnd < ss "
			if avg == 0:
			    avg = 'NA'
		if avg != 'NA':
		    mcollection.update({"_id":test["_id"]},{"$set":{"cons."+consName: avg}}, upsert = True)
    return 0

def updateMap(path, db, motifChrom="chr17", window=100):
    """update mapability scores from bedGraph files"""
    mcollection = db["hg19"+motifChrom]
    for infile in glob.glob(os.path.join(path, "wgEncodeCrgMapabilityAlign24mer"+motifChrom+".bedGraph.bz2")):##better in bedGraph
	(mappath,mapfilename) = os.path.split(infile)
	expName = mapfilename.split(motifChrom)[0]
	print 'updating', expName, motifChrom
	#print expName
	with bz2.BZ2File(infile,'r') as bedGraphFile:
	    #bedGraph = csv.reader(bedGraphFile, delimiter = '\t')
	    bbFile = os.path.join(mappath,expName+motifChrom+'.bb')
	    if not os.path.isfile(bbFile):
		countBed.compressBed4(bedGraphFile, expName, bbFile)
	    coordDict, valuesDict = countBed.getBinBedCoord(bbFile, expName)
	    arrayDict=defaultdict(list)
	    cursor = mcollection.find({"tf_name":{"$in": ["IRF3","MAFK","NFYA","SIN3A","ZNF384"]}})
		#{"tf_name":tfName})
	    for test in cursor:
		motifStart, motifEnd = test["genomic_region"]["start"], \
                                       test["genomic_region"]["end"]
		if not motifChrom in arrayDict:
		    arrayDict[motifChrom] = countBed.buildBedHist(motifChrom, coordDict, valuesDict, expName)
		xs, xvals, sums = arrayDict[motifChrom]
		avg = countBed.queryHist(xs, xvals, sums, motifStart-window, motifEnd+window)[0]		
		#if avg != 0:
		    #print motifChrom, motifStart, motifEnd
		    #print avg
 		mcollection.update({"_id":test["_id"]},{"$set":{"map."+expName:avg}}, upsert = True)
    return 0

def updateChip(path,db,motifChrom="chr17",window=0):
    """mark motifs if it overlaps ENCODE excluded regions"""
    mcollection = db["hg19"+motifChrom]
    for infile in glob.glob(os.path.join(path,"*.bed.gz")):
	(regpath,regfilename) = os.path.split(infile)
	expName = regfilename.split(motifChrom)[0]
	with gzip.open(infile,'rt') as bedFile:
	    bed = csv.reader(bedFile, delimiter = '\t') 
	    annoIntvlDict = countBed.getBed4Anno(bed,expName)
	    intervalStartDict = countBed.sortStart(annoIntvlDict)
	    intervalEndDict = countBed.sortEnd(annoIntvlDict)
	    cursor = mcollection.find({"tf_name":{"$in": ["IRF3","MAFK","NFYA","SIN3A","ZNF384"]}})
		#{"tf_name": tfName})
	    for test in cursor:
		motifStart, motifEnd = test["genomic_region"]["start"], \
					test["genomic_region"]["end"]
	        regionList, valueList = countBed.getMotifAnno(annoIntvlDict,
				intervalStartDict,intervalEndDict,motifChrom,
				motifStart,motifEnd,window)
		if valueList != []:
			#print regionList, valueList, motifChrom, motifStart, motifEnd
			mcollection.update({"_id": test["_id"]}, 
			{"$set": {"chip."+expName: valueList[0][1]}}, upsert = True)
	    	
    return 0

def updateOvlpRegions(path,db,motifChrom='chr17',window=0, stranded=True):
    mcollection = db["hg19"+motifChrom]
    for infile in glob.glob(os.path.join(path,"*.bed.gz")):
	expName="exon1st"
	#(regpath,regfilename) = os.path.split(infile)
	#expName = regfilename.split('.')[0]
	print 'updating', expNam
	with gzip.open(infile,'rt') as bedFile:
	    bed = csv.reader(bedFile, delimiter='\t')
	    annoIntvlDict = countBed.getBed6Anno(bed,expName)
	    intervalStartDict = countBed.sortStart(annoIntvlDict)
	    intervalEndDict = countBed.sortEnd(annoIntvlDict)
	    cursor = mcollection.find({"tf_name":{"$in": ["IRF3","MAFK","NFYA","SIN3A","ZNF384"]}})
	    if stranded:
	        for test in cursor:
		    exonList = []
		    motifStart, motifEnd, motifStrand = test["genomic_region"]["start"], \
					test["genomic_region"]["end"], \
					test["genomic_region"]["strand"]
		    regionList, valueList = countBed.getMotifAnno(annoIntvlDict,
				intervalStartDict,intervalEndDict,motifChrom,
				motifStart,motifEnd,window)
		    if valueList != []:##[(exon_name,strand),]
			for i in xrange(len(valueList)):
			    if valueList[i][1] == motifStrand:
				exonList.append(regionList[i])
		    if exonList != []:
			mcollection.update({"_id": test["_id"]},
				{"$set":{expName: exonList}}, upsert = True)
	    else:
		for test in cursor:
		    motifStart, motifEnd = test["genomic_region"]["start"], \
				test["genomic_region"]["end"]
		    regionList, valueList = countBed.getMotifAnno(annoIntvlDict,
				intervalStartDict,intervalEndDict,motifChrom,
				motifStart,motifEnd,window)
		    if regionList != []:
			mcollection.update({"_id": test["_id"]},
				{"$set":{"map."+expName: regionList}}, upsert = True)
    return 0


def updateMotifGenomicRegions(path,db,chromNum,mcollection=None):
    """clean overlapping motif entries and update motif genomic loci"""
    for infile in glob.glob(os.path.join(path,"hg19_"+chromNum+"_1e3_100_*_fimoout_092914_sort.txt.gz")):
        (fimopath,fimofilename) = os.path.split(infile)
        tfName = fimofilename.split('hg19_'+chromNum+'_1e3_100_')[-1].split('fimoout')[0].split('_')[0]
	print "updating", tfName
        f = gzip.open(infile)
        fimo = csv.reader(f,delimiter='\t')
        pchrom = 'chr0'
        rowlist = []
        for row in fimo:
            if not "pattern" in row[0]:
                ##fimo output is 1-based, so don't have to change anything in mongodb or gff
                chrom, start, end, strand, pval = row[1], int(row[2]), int(row[3]), row[4], float(row[6])
		if mcollection is None:
		    print 'initiating collection', chrom
		    mcollection = db["hg19"+chrom]
		    mcollection.ensure_index("tf_name",name="tf_name",unique=False,background=True)
		    mcollection.ensure_index("motif_id",name="motif_id",unique=False,background=True)
		    mcollection.ensure_index([
				("genomic_region.start",DESCENDING),
				("genomic_region.end",ASCENDING),
				("genomic_region.strand",DESCENDING)],
				name="genomic_region",unique=False,background=True)
		    mcollection.ensure_index("motif_score",name="motif_score",unique=False,background=True)
		    print 'initial collection count', mcollection.count()
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
                            #m = mcollection.find_one({"motif_id": wmid, "tf_name": tfName})
                            mcollection.insert({"motif_id": wmid, 
						"tf_name": tfName, 
						#"motif_gene_mapping_info": m["motif_gene_mapping_info"],
						#"motif_ct_info": m["motif_ct_info"], 
						#"motif_tf_info": m["motif_tf_info"],
						#"motif_mapability_info": m["motif_mapability_info"],
						#"motif_cons_info": m["motif_cons_info"], 
						"genomic_region": 
						{"chr": wchrom, 
						"start": wstart, 
						"end": wend, 
						"strand": wstrand}, 
						"motif_score": wpval})
                            #mcollection.remove({"tf_name": tfName, "motif_score": None})##remove original placeholder entry
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
                        #m = mcollection.find_one({"motif_id": wmid, "tf_name": tfName})
                        mcollection.insert({"motif_id": wmid, 
					    "tf_name": tfName, 
					    #"motif_gene_mapping_info": m["motif_gene_mapping_info"],
					    #"motif_ct_info": m["motif_ct_info"], 
					    #"motif_tf_info": m["motif_tf_info"], 
					    #"motif_mapability_info": m["motif_mapability_info"],
					    #"motif_cons_info": m["motif_cons_info"],
					    "genomic_region": 
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
    if len(argv) < 4:
        sys.stderr.write("Usage: %s motif_tf_info_file path-to-fimo-output chromNum\n" % argv[0])
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

    #global mcollection 
    #mcollection = db["hg19chr17"]
    ##drop collection
    #c["mm9"].drop_collection("motif_instance_hughes_test")
    #mcollection.remove()
    #db.drop_collection('hg19chr17')
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
    chromNum = sys.argv[3]

    ifile = open(infile,'rt')
    tf_info = csv.reader(ifile, delimiter = '\t')

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

    #cursor = mcollection.find()
    #print 'before entering genomic region info ', mcollection.count()
    #c = iter(cursor)

    if chromNum[-2] == '0':
	motifChrom = 'chr'+chromNum[-1]
    elif not chromNum[5:] == '20':
	motifChrom = 'chr'+chromNum[5:]
    else:
	motifChrom = 'chrX'

    db = client["hg19"+motifChrom]
    ##clean overlapping motif entries
    startTime = time.clock()
    #updateMotifGenomicRegions(path, db, chromNum)
    #updateCount(path,db,motifChrom)
    #updateGC(path,db,motifChrom)
    #updateExcludedRegions(path,db)
    #updatePi(path,db,motifChrom)
    #updateMap(path,db,motifChrom)
    #updateCons(path,db,motifChrom)

    ##write gff
#   print 'updated count ', cursor.count(), 'total count after update', mcollection.count(), 'update time', time.clock() - startTime

    updateFOS(path, db, motifChrom)#, stranded=True)#'Ctcf', 'chr17', method="Binom")
    #ofile = open('/home/xc406/data/mongodbtest/ctcfData.txt','wt')
    #writer = csv.writer(ofile, delimiter='\t')

    ##update gene features
#    refSeqFile = open('/home/xc406/data/hg19_refseq_June_2014.txt','rt')
#    refSeqReader = csv.reader(refSeqFile, delimiter='\t')
#    tssDict, geneNameDict, geneRangeDict = getRefSeqDict(refSeqReader)
    #mcollection = db["hg19"+motifChrom]
    #cursor = mcollection.find()#{"tf_name":"Stat3"})
    #intervalDict = sortInterval(geneRangeDict)

#    intervalStartDict = countBed.sortStart(geneRangeDict)
#    intervalEndDict = countBed.sortEnd(geneRangeDict)
    #with client.start_request():##open update
    #print "Zic1 motif: {0}".format(test)
#    motifChrom = 'chr17'
#    for test in cursor:
#        motifStart, motifEnd = test["genomic_region"]["start"], test["genomic_region"]["end"]
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
#	mcollection.update({"_id":test["_id"]},{"$set":{"target_gene":{"closest": (closest[0],closest[2]), 
#			"dist_tss": closest[1], "10kb": t[1]}}}, upsert = True)
    	#test["genomic_regions_gene_mapping"]["closest_gene"] = (closest[0],closest[2])#closestGene(tssDict,geneNameDict,motifChrom,motifStart,motifEnd)[0][0] 
   	#test["genomic_regions_gene_mapping"]["dist_tss"] = closest[1]#closestGene(tssDict,geneNameDict,motifChrom,motifStart,motifEnd)[1]
	#test["genomic_regions_gene_mapping"]["genelist10kb"] = t[0]
	#test["genomic_regions_gene_mapping"]["transcriptidlist10kb"] = t[1]#getTargetGene(geneRangeDict,geneNameDict,motifChrom, motifStart, motifEnd, 0)
    	#mcollection.save(test)
	#print "Hes5 motif: {0}".format(test)	

    #print 'update time', time.time() - startTime2

    #writeMongo.getData(cursor,writer,0)
    #ofile.close()

    #print "Zic1 motif update: {0}".format(testupdate)
    print 'total time', time.clock() - startTime
    #mcollection.insert({motif_id: item['motif_id']})
    #c = cursor.next()
    #print cursor.next()
    ##update motif_genomic_regions_info
    #db.motif_instance_hughes_test.insert({"motif_tf_info": db.motif_instance_hughes_test.findOne({'tf_name':'Irf4'}).motif_tf_info})

if __name__=='__main__':
    sys.exit(main(sys.argv))
