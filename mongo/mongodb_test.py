from pymongo import MongoClient
from pymongo.errors import DuplicateKeyError
import os, csv, sys, glob
from pymongo import ASCENDING, DESCENDING
from pymongo import Connection
import numpy as np
from bisect import bisect
from collections import defaultdict

def getTssDict(refSeqReader):
    """store transcripts loci info in dictionaries"""
    tssDict = defaultdict(list)
    geneNameDict = defaultdict(str)
    sizeDict = defaultdict(int)
    for row in refSeqReader:
	try:
	    chrom, strand, txStart, txEnd, geneName = row[2], row[3], int(row[4]), int(row[5]), row[12]
	    size = txEnd - txStart
	    if strand == '+':
	        if not geneName in sizeDict:
	    	    tssDict[chrom].append(txStart)
	    	    geneNameDict[(chrom,txStart)] = geneName
		    sizeDict[geneName] = size
	        else:
		    if size > sizeDict[geneName]:
		    	tssDict[chrom].append(txStart)
		    	geneNameDict[(chrom,txStart)] = geneName
		    	sizeDict[geneName] = size
	    else:
	    	if not geneName in sizeDict:
		    tssDict[chrom].append(txEnd)
		    geneNameDict[(chrom,txEnd)] = geneName
		    sizeDict[geneName] = size
	    	else:
		    if size > sizeDict[geneName]:
	    	    	tssDict[chrom].append(txEnd)
	    	    	geneNameDict[(chrom,txEnd)] = geneName
		    	sizeDict[geneName] = size
	except ValueError:
	    pass
    return tssDict, geneNameDict

def closestGene(tssDict,geneNameDict,motifChrom,motifStart,motifEnd):
    """find closest gene feature"""
    for chrom in tssDict:
    	tssDict[chrom].sort()##don't sort beforehand, change the original dictionary
    motifMid = motifStart+(motifEnd-motifStart)/2
    i = bisect(tssDict[motifChrom],motifMid)
    if i < len(tssDict[motifChrom]) and i > 0:
    	if abs(motifMid - tssDict[motifChrom][i-1]) >= abs(tssDict[motifChrom][i] - motifMid):
    	    geneName = geneNameDict[(motifChrom,tssDict[motifChrom][i])]
	    dist = abs(tssDict[motifChrom][i] - motifMid)
    	else:
	    geneName = geneNameDict[(motifChrom,tssDict[motifChrom][i-1])]
	    dist = abs(motifMid - tssDict[motifChrom][i-1])
    elif i == len(tssDict[motifChrom]):
	if not i == 0:
	    geneName = geneNameDict[(motifChrom,tssDict[motifChrom][i-1])]
	    dist = abs(motifMid - tssDict[motifChrom][i-1])
	else:
	    geneName = None
	    dist = None
    elif i == 0:
	geneName = geneNameDict[(motifChrom,tssDict[motifChrom][i])]
	dist = abs(motifMid - tssDict[motifChrom][i])
    return geneName, dist

def makeGff(cursor, gffWriter, window):
    """write gff files with a mongodb query"""
    for item in cursor:
	#print item["genomic_regions_gene_mapping"]
	info = 'gene_id 1' #+ item["genomic_regions_gene_mapping"]["closest_gene"] + '; ct_name ' + item["ct_info"]["ct_name"] #+ '; scores ' + "_".join(item["ct_info.accessibility_score"])
	row = [item["motif_genomic_regions_info"]["chr"],item["tf_name"],"motif",
		item["motif_genomic_regions_info"]["start"]-window,
		item["motif_genomic_regions_info"]["end"]+window,
		item["motif_score"],item["motif_genomic_regions_info"]["strand"],'.',info]
	gffWriter.writerows([row])

def getData(mcollection, ofile, tfName):
    cursor = mcollection.find({"tf_name": tfName})
    for m in cursor:
	row = [m["genomic_region"]["chr"],
		m["genomic_region"]["start"],
		m["genomic_region"]["end"],
		m["motif_score"],
		m["target_gene"]["dist_tss"],
		m["cons"]["pi"],
		m["cons"]["piSnp"],
		m["cons"]["piIndel"],
		m["cons"]["phastCons46way_primates"],
		m["cons"]["phyloP46way_primate"],
		m["cons"]["phastCons100way_wigFix"],
		m["cons"]["phyloP100way_wigFix"],
		m["gc"][0],
		m["gc"][1],
		m["mapability"]["wgEncodeCrgMapabilityAlign24mer"]]
	writer.writerows([row])


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
    db = client["test"]
    mcollection = db["chr17"]
    ##drop collection
    #mcollection.remove()
    #db.drop_collection('motif_instance_hughes_test')
    #print 'clean', mcollection.count()
    #mcollection.ensure_index("motif_id",name="m_id",unique=False,background=True)
    #mcollection.ensure_index("tf_name",name="tf_name",unique=False,background=True)
    #collection.ensure_index("motif_type",name="motif_type",unique=True,drop_dups=True,background=True)
    #index genomic regions
    #mcollection.ensure_index([("motif_genomic_regions_info.chr",DESCENDING),
#				("motif_genomic_regions_info.start",DESCENDING),
#				("motif_genomic_regions_info.end",ASCENDING),
#				("motif_genomic_regions_info.strand",DESCENDING)],
#				name="genomic_reg",unique=False,background=True)
    #collection.ensure_index("motif_id",name="motif_id",unique=True,drop_dups=True,background=True)
	
    infile = sys.argv[1]#'/home/xc406/data/hg19motifs90/TF_Information90hg19.txt'
    path = sys.argv[2]

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
		"genomic_regions_gene_mapping":{
			"closest_gene": None,
			"feature": None,#intergenic or 3' 5'
			"dist_tss": None,
			"genelist10kb": []##list of gene--center of motif fall in gene plus and minus 10kb
			#"epu_id": #boolean		
		},
		"ct_info":{
			"ct_name": None,
			"ct_type": None, #normal/cancerous/cellline/primarycell
			"accessibility_score": {}, #type: log likelihood score
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

#        try:	
		#print motif_instance	
#		mcollection.insert(motif_instance)
#		del motif_instance
		#print 'inserted' #collection
#    	except DuplicateKeyError:
		#print 'dup'
#		pass

    #cursor = mcollection.find()
    #print 'before entering genomic region info ', mcollection.count()
    #c = iter(cursor)

    ##clean overlapping motif entries
    #updateMotifLoci(path)

    ##write gff
#    cursor = mcollection.find({"tf_name": "Hes5"})
#    print 'updated count ', cursor.count()#mcollection.count()
#    ofile = open('/home/xc406/data/mongodbtest/test.gff','wt')
#    gffWriter = csv.writer(ofile, delimiter='\t')
#    makeGff(cursor,gffWriter,0)

    ##update gene features
    refSeqFile = open('/home/xc406/data/mm9_refseq_June_2014.txt','rt')
    refSeqReader = csv.reader(refSeqFile, delimiter='\t')
    tssDict, geneNameDict = getTssDict(refSeqReader)
    cursor = mcollection.find({"tf_name":"Hes5"})
    #print "Zic1 motif: {0}".format(test)
    for test in cursor:
        motifChrom, motifStart, motifEnd = test["motif_genomic_regions_info"]["chr"], test["motif_genomic_regions_info"]["start"], test["motif_genomic_regions_info"]["end"]
    	#print closestGene(tssDict,geneNameDict,motifChrom,motifStart,motifEnd)[0][0]
    	#print closestGene(tssDict,geneNameDict,motifChrom,motifStart,motifEnd)[1]
    	#print test["motif_id"], test["motif_genomic_regions_info"]["chr"],test["motif_genomic_regions_info"]["start"]
#    mcollection.update({"motif_id": test["motif_id"], "motif_genomic_regions_info":{"chr": test["motif_genomic_regions_info"]["chr"], 
#			"start": test["motif_genomic_regions_info"]["start"]}},{"$set": 
    #mcollection.update(test,{"$set":{"genomic_regions_gene_mapping":{"closest_gene": closestGene(tssDict,geneNameDict,motifChrom,motifStart,motifEnd)[0][0], 
#			"dist_tss": closestGene(tssDict,geneNameDict,motifChrom,motifStart,motifEnd)[1]}}})
	#if len(closestGene(tssDict,geneNameDict,motifChrom,motifStart,motifEnd)[0]) > 1:
	    #print closestGene(tssDict,geneNameDict,motifChrom,motifStart,motifEnd)[0], motifChrom, motifStart, motifEnd	    
    	test["genomic_regions_gene_mapping"]["closest_gene"] = closestGene(tssDict,geneNameDict,motifChrom,motifStart,motifEnd)[0][0] 
    	test["genomic_regions_gene_mapping"]["dist_tss"] = closestGene(tssDict,geneNameDict,motifChrom,motifStart,motifEnd)[1]
    	mcollection.save(test)
	
    #print cursor.next()

if __name__=='__main__':
    sys.exit(main(sys.argv))
