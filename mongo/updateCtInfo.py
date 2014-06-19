from pymongo import MongoClient
from pymongo.errors import DuplicateKeyError
import os, csv, sys, glob
from pymongo import ASCENDING, DESCENDING
from pymongo import Connection
import countWig
from collections import defaultdict

def getCount(path, tfName):
    for infile in glob.glob(os.path.join(path,"*.wig")):
        (wigpath,wigfilename) = os.path.split(infile)
        ##depends on the data type and source
	methodName = "Dnase"
        ctName = wigfilename.split('EncodeUwDnase')[-1].split('Aln')[0]
        wigFile = open(infile,'rt')
        wig = csv.reader(wigFile,delimiter='\t')
        coordDict, valuesDict = countWig.getCoord(wig,ctName)
        arrayDict = defaultdict(list)
        cursor = mcollection.find({"tf_name": tfName})
	for test in cursor:
            motifChrom, motifStart, motifEnd = test["motif_genomic_regions_info"]["chr"], test["motif_genomic_regions_info"]["start"], test["motif_genomic_regions_info"]["end"]
   	    if not motifChrom in arrayDict:
	    	arrayDict[motifChrom] = countWig.buildHist(motifChrom,coordDict,valuesDict,ctName)
            xs, xvals, sums = arrayDict[motifChrom]
            count = countWig.queryHist(xs, xvals, sums, motifStart, motifEnd)[0]
	    print count
	    test["ct_info"]["accessibility_score"][methodName] = count
	    mcollection.save(test)
    return 0

def main(argv):
    if len(argv) < 3:
        sys.stderr.write("Usage: %s infile path-to-wig-files\n" % argv[0])
        return 1
    if not os.path.isfile(argv[1]):
        sys.stderr.write('Error: infile %r was not found!\n' % argv[1])
        return 1
    if not os.path.exists(argv[2]):
        sys.stderr.write('Error: path-to-wig-files %r was not found!\n' % argv[2])
        return 1 

    server = 'localhost'
    port = 27017
    client = MongoClient(server, port)
    db = client["mm9"]
    global mcollection
    mcollection = db["motif_instance_hughes_test"]
    ##drop collection
    #c = Connection()
    #c["mm9"].drop_collection("motif_instance_hughes_test")
    #mcollection.remove()
    #db.drop_collection('motif_instance_hughes_test')
    #mcollection = db["motif_instance_hughes_test"]
    #print 'clean', mcollection.count()
    #mcollection.ensure_index("motif_id",name="m_id",unique=False,background=True)
    #mcollection.ensure_index("tf_name",name="tf_name",unique=False,background=True)
    #collection.ensure_index("motif_type",name="motif_type",unique=True,drop_dups=True,background=True)
    #index genomic regions
    #mcollection.ensure_index([("motif_genomic_regions_info.chr",DESCENDING),("motif_genomic_regions_info.start",DESCENDING),("motif_genomic_regions_info.end",ASCENDING),("motif_genomic_regions_info.strand",DESCENDING)],name="genomic_reg",unique=False,background=True)
    #collection.ensure_index("motif_id",name="motif_id",unique=True,drop_dups=True,background=True)
    #print collection
	
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
    #    try:	
		#print motif_instance	
#		mcollection.insert(motif_instance)
#		del motif_instance
		#print 'inserted' #collection
 #   	except DuplicateKeyError:
		#print 'dup'
#		pass
 #   cursor = mcollection.find()
    print 'before entering genomic region info ', mcollection.count()
#    c = iter(cursor)
    ##clean overlapping motif entries
    getCount(path,"Hes5")
#    cursor = mcollection.find()
    print 'updated count ', mcollection.count()
    #for item in cursor:
	#mcollection.insert({motif_id: item['motif_id']})
	#print item
    #    print item['tf_name']
	#c = cursor.next()
	#print cursor.next()	
	#collection.find_one({"motif_id":}) 





if __name__=='__main__':
    sys.exit(main(sys.argv))
