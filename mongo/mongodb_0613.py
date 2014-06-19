from pymongo import MongoClient
from pymongo.errors import DuplicateKeyError
import os, csv, sys, glob
from pymongo import ASCENDING, DESCENDING
from pymongo import Connection
import numpy as np
from bisect import bisect
from collections import defaultdict
import time
from ngs_plumbing.intervals import Interval, IntervalList

def updateMotifGenomicRegions(path):
    ##clean overlapping motif entries
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
                ##########strand specific############ if control for strand-specificity, split fimo files as + and -, and run separately
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
                            #mcollection.update({"motif_id": wmid, "tf_name": tfName},
				#{"$set": {"motif_genomic_regions_info.chr": wchrom, 
				#"motif_genomic_regions_info.start": start, "motif_genomic_regions_info.end": end, 
				#"motif_genomic_regions_info.strand": strand, 
				#"motif_score": wpval, "genomic_regions_gene_mapping": "closest_gene": None, 
				#"feature": None, "dist_tss": None} }, upsert = True)
                            #print "Found m: {0}".format(m)
                            mcollection.insert({"motif_id": wmid, 
						"tf_name": tfName, 
						"motif_gene_mapping_info": m["motif_gene_mapping_info"],
						"motif_ct_info": m["motif_ct_info"], 
						"motif_tf_info": m["motif_tf_info"],
						"motif_mapability_info": m["motif_mapability_info"],
						"motif_cons_info": m["motif_cons_info"], 
						"motif_genomic_regions_info": 
						{"chr": wchrom, 
						"start": wstart, 
						"end": wend, 
						"strand": wstrand}, 
						"motif_score": wpval})
                            mcollection.remove({"tf_name": tfName, "motif_score": None})##remove original placeholder entry
                            #print mcollection.count()
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
					    "motif_gene_mapping_info": m["motif_gene_mapping_info"],
					    "motif_ct_info": m["motif_ct_info"], 
					    "motif_tf_info": m["motif_tf_info"], 
					    "motif_mapability_info": m["motif_mapability_info"],
					    "motif_cons_info": m["motif_cons_info"],
					    "motif_genomic_regions_info": 
					     {"chr": wchrom, 
					      "start": wstart, 
					      "end": wend, 
					      "strand": wstrand}, 
					    "motif_score": wpval})
                        #print mcollection.count()
                        #mcollection.update({"motif_id": wmid, "tf_name": tfName}, 
				#{"$set": {"motif_genomic_regions_info.chr": wchrom, 
				#"motif_genomic_regions_info.start": start, "motif_genomic_regions_info.end": end, 
				#"motif_genomic_regions_info.strand": strand, "motif_score": wpval}})
                        #ofile.writelines('\t'.join(modRow(wrow))+'\n')
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
		    geneRangeDict[chrom][Interval(txStart,txEnd)] = (geneName, transcriptId)
	        else:
		    if size > sizeDict[geneName]:
		    	tssDict[chrom].append(txStart)
		    	geneNameDict[(chrom,txStart)] = (geneName, transcriptId)
		    	sizeDict[geneName] = size
			geneRangeDict[chrom][Interval(txStart,txEnd)] = (geneName, transcriptId)
	    else:
	    	if not geneName in sizeDict:
		    tssDict[chrom].append(txEnd)
		    geneNameDict[(chrom,txEnd)] = (geneName, transcriptId)
		    sizeDict[geneName] = size
		    geneRangeDict[chrom][Interval(txStart,txEnd)] = (geneName, transcriptId)
	    	else:
		    if size > sizeDict[geneName]:
	    	    	tssDict[chrom].append(txEnd)
	    	    	geneNameDict[(chrom,txEnd)] = (geneName, transcriptId)
		    	sizeDict[geneName] = size
			geneRangeDict[chrom][Interval(txStart,txEnd)] = (geneName, transcriptId)
	except ValueError:
	    pass
    return tssDict, geneNameDict, geneRangeDict

def closestGene(tssDict,geneNameDict,motifChrom,motifStart,motifEnd):
    """find closest gene feature"""
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

def sortInterval(geneRangeDict):
    """sort all intervals on a chromosome"""
    intervalDict = {}
    for chrom in geneRangeDict:
	itl = IntervalList( intv for intv in geneRangeDict[chrom])
	itl.sort()
	intervalDict[chrom] = itl
    return intervalDict
        
def getTargetGene(geneRangeDict,intervalDict,motifChrom, motifStart, motifEnd, txWindow):
    """push motif target genes into a list"""
    targetList = []
    transcriptList = []
    #itl = IntervalList( intv for intv in geneRangeDict[motifChrom])
    #itl.sort()
    try:
    	intervals = intervalDict[motifChrom]
    	motifInterval = Interval(motifStart, motifEnd)
    	overlapping = [ x for x in intervals if (x[1]>motifInterval[0] and x[0]<motifInterval[0]) or (x[1]>motifInterval[1] and x[0]<motifInterval[1]) ]
        #print 'overlapping',overlapping
    	if not len(overlapping) == 0:
	    for x in overlapping:
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

def makeGff(cursor, gffWriter, window):
    """write gff files with a mongodb query"""
    for item in cursor:
	#print item["genomic_regions_gene_mapping"]
	info = 'gene_id ' + '_'.join(item["motif_gene_mapping_info"]["genelist10kb"])#+ 
		#'; ct_name ' + item["ct_info"]["ct_name"] #+ '; scores ' + "_".join(item["ct_info.accessibility_score"])
	row = [item["motif_genomic_regions_info"]["chr"],
		item["tf_name"],
		"motif",
		item["motif_genomic_regions_info"]["start"]-window,
		item["motif_genomic_regions_info"]["end"]+window,
		item["motif_score"],
		item["motif_genomic_regions_info"]["strand"],
		'.',
		info]
	gffWriter.writerows([row])
	
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
        sys.stderr.write("Usage: %s infile path-to-fimo-output\n" % argv[0])
        return 1
    if not os.path.isfile(argv[1]):
        sys.stderr.write('Error: infile %r was not found!\n' % argv[1])
        return 1
    if not os.path.exists(argv[2]):
        sys.stderr.write('Error: path-to-fimo-output %r was not found!\n' % argv[2])
        return 1 

    server = 'localhost'
    port = 27017
    client = MongoClient(server, port)
    db = client["mm9"]
    global mcollection 
    mcollection = db["motif_instance_hughes_90"]
    ##drop collection
    #c = Connection()
    #c["mm9"].drop_collection("motif_instance_hughes_test")
    #mcollection.remove()
    #db.drop_collection('motif_instance_hughes_test')
    #mcollection = db["motif_instance_hughes_test"]
    #print 'clean', mcollection.count()
#    mcollection.ensure_index("motif_id",name="m_id",unique=False,background=True)
#    mcollection.ensure_index("tf_name",name="tf_name",unique=False,background=True)
    #collection.ensure_index("motif_type",name="motif_type",unique=True,drop_dups=True,background=True)
#    mcollection.ensure_index("motif_gene_mapping_info.genelist10kb",name = "target_gene",unique=False,background=True)
#    mcollection.ensure_index("motif_gene_mapping_info.closest_gene",name = "closest_gene",unique=False,background=True)
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
		"motif_mapability_info":{
				"exclude": [],
				"score": [],
				"gc_content": None
		},
		"motif_cons_info":{
				"phylop_primates": None,
				"phylop_mammals": None,
				"phylop_vertebrate": None,
				"phastCons_primates": None,
				"phastCons_mammals": None,
				"phastCons_vertebrate": None
		},
		"motif_gene_mapping_info":{
			"closest_gene": None,
			"feature": None,#intergenic or 3' 5'
			"dist_tss": None,
			"genelist10kb": [],##list of gene--center of motif fall in gene plus and minus 10kb
			"transcriptidlist10kb": []
			#"epu_id": #boolean		
		},
		"motif_ct_info":{
			"ct_name": None,
			"ct_type": None, #normal/cancerous/cellline/primarycell
			"accessibility_score": {}, #{type: log likelihood score}
			#"accessibility_type": [], #dhs, dgf, faire
			"chip_score": {}, #tf: pval(overlapping peaks)
			"h3k4me3_score": None,
			"h3k4me1_score": None,
			"h3k27ac_score": None,
			"p300_score": None,
			"pol2_score": None
		},	
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
    print 'before entering genomic region info ', mcollection.count()
    #c = iter(cursor)

    ##clean overlapping motif entries
    startTime = time.time()
    #updateMotifGenomicRegions(path)

    ##write gff
    #cursor = mcollection.find({"tf_name": "Stat3"})
    #print 'updated count ', cursor.count(), 'total count after update', mcollection.count(), 'update time', time.time() - startTime
    #ofile = open('/home/xc406/data/mongodbtest/test.gff','wt')
    #gffWriter = csv.writer(ofile, delimiter='\t')
#    makeGff(cursor,gffWriter,0)

    ##update gene features
    refSeqFile = open('/home/xc406/data/mm9_refseq_June_2014.txt','rt')
    refSeqReader = csv.reader(refSeqFile, delimiter='\t')
    tssDict, geneNameDict, geneRangeDict = getRefSeqDict(refSeqReader)
    cursor = mcollection.find({"tf_name":"Stat3"})
    intervalDict = sortInterval(geneRangeDict)
    startTime2 = time.time()
    #print "Zic1 motif: {0}".format(test)
    for test in cursor:
        motifChrom, motifStart, motifEnd = test["motif_genomic_regions_info"]["chr"], test["motif_genomic_regions_info"]["start"], test["motif_genomic_regions_info"]["end"]
    	#print closestGene(tssDict,geneNameDict,motifChrom,motifStart,motifEnd)[0][0]
    	#print closestGene(tssDict,geneNameDict,motifChrom,motifStart,motifEnd)[1]
    	#print test["motif_id"], test["motif_genomic_regions_info"]["chr"],test["motif_genomic_regions_info"]["start"]
	#startTime = time.time()
	t = getTargetGene(geneRangeDict,intervalDict,motifChrom, motifStart, motifEnd, 10000)
	#endTime1 = time.time()
	#print "all target mapping", t, motifChrom, motifStart, motifEnd, endTime1 - startTime
	closest = closestGene(tssDict,geneNameDict,motifChrom,motifStart,motifEnd)
	#endTime2 = time.time()
	#print "closest gene mapping time", closest, endTime2 - endTime1 
#    mcollection.update({"motif_id": test["motif_id"], "motif_genomic_regions_info":{"chr": test["motif_genomic_regions_info"]["chr"], 
#			"start": test["motif_genomic_regions_info"]["start"]}},{"$set": 
    #mcollection.update(test,{"$set":{"genomic_regions_gene_mapping":{"closest_gene": closestGene(tssDict,geneNameDict,motifChrom,motifStart,motifEnd)[0][0], 
#			"dist_tss": closestGene(tssDict,geneNameDict,motifChrom,motifStart,motifEnd)[1]}}})
	#if len(closestGene(tssDict,geneNameDict,motifChrom,motifStart,motifEnd)[0]) > 1:
	    #print closestGene(tssDict,geneNameDict,motifChrom,motifStart,motifEnd)[0], motifChrom, motifStart, motifEnd
    	test["motif_gene_mapping_info"]["closest_gene"] = (closest[0],closest[2])#closestGene(tssDict,geneNameDict,motifChrom,motifStart,motifEnd)[0][0] 
   	test["motif_gene_mapping_info"]["dist_tss"] = closest[1]#closestGene(tssDict,geneNameDict,motifChrom,motifStart,motifEnd)[1]
	test["motif_gene_mapping_info"]["genelist10kb"] = t[0]
	test["motif_gene_mapping_info"]["transcriptidlist10kb"] = t[1]#getTargetGene(geneRangeDict,geneNameDict,motifChrom, motifStart, motifEnd, 0)
    	mcollection.save(test)
	#print "Hes5 motif: {0}".format(test)	

    print 'update time', time.time() - startTime2
    #cursor = mcollection.find({"tf_name":"Zic1"})
    #makeGff(cursor,gffWriter,0)
    #ofile.close()
    #testupdate = mcollection.find_one({"tf_name":"Zic1"})#{"motif_id": test["motif_id"], 
		#"motif_genomic_regions_info":{"chr": test["motif_genomic_regions_info"]["chr"], "start": test["motif_genomic_regions_info"]["start"]}})
    #print "Zic1 motif: {0}".format(test)
    #print "Zic1 motif update: {0}".format(testupdate)
    print 'total time', time.time() - startTime
    #mcollection.insert({motif_id: item['motif_id']})
    #c = cursor.next()
    #print cursor.next()
    ##update motif_genomic_regions_info
    #db.motif_instance_hughes_test.insert({"motif_tf_info": db.motif_instance_hughes_test.findOne({'tf_name':'Irf4'}).motif_tf_info})

if __name__=='__main__':
    sys.exit(main(sys.argv))
