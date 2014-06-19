from pymongo import MongoClient
from pymongo.errors import DuplicateKeyError
import os, csv, sys, glob
from pymongo import ASCENDING, DESCENDING
from pymongo import Connection

def modRow(row):
    """fix motif size in gff file"""
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
    mcollection = db["motif_instance_hughes_test"]
    ##drop collection
    #c = Connection()
    #c["mm9"].drop_collection("motif_instance_hughes_test")
    #mcollection.remove()
    db.drop_collection('motif_instance_hughes_test')
    mcollection = db["motif_instance_hughes_test"]
    print 'clean', mcollection.count()
    mcollection.ensure_index("motif_id",name="m_id",unique=False,background=True)
    mcollection.ensure_index("tf_name",name="tf_name",unique=False,background=True)
    #collection.ensure_index("motif_type",name="motif_type",unique=True,drop_dups=True,background=True)
    #index genomic regions
    mcollection.ensure_index([("motif_genomic_regions_info.chr",DESCENDING),("motif_genomic_regions_info.start",DESCENDING),("motif_genomic_regions_info.end",ASCENDING),("motif_genomic_regions_info.strand",DESCENDING)],name="genomic_reg",unique=False,background=True)
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
        try:	
		#print motif_instance	
		mcollection.insert(motif_instance)
		del motif_instance
		#print 'inserted' #collection
    	except DuplicateKeyError:
		#print 'dup'
		pass
    cursor = mcollection.find()
    print 'before entering genomic region info ', mcollection.count()
    c = iter(cursor)
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
				##best of the overlapping motifs, current line slightly overlapping but ok
			    elif len(rowlist) == 1:
			        wrow = prow
				##previous line non-overlapping, current line slightly overlapping but ok
			    rowlist = []
			    rowlist.append(row)
			    wmid, wchrom, wstart, wend, wstrand, wpval = wrow[0], wrow[1], int(wrow[2]), int(wrow[3]), wrow[4], float(wrow[6])
			    m = mcollection.find_one({"motif_id": wmid, "tf_name": tfName})
			    #if m is None:
				#mcollection.update({"motif_id": wmid, "tf_name": tfName},{"$set": {"motif_genomic_regions_info.chr": wchrom, "motif_genomic_regions_info.start": start, "motif_genomic_regions_info.end": end, "motif_genomic_regions_info.strand": strand, "motif_score": wpval, "genomic_regions_gene_mapping": "closest_gene": None, "feature": None, "dist_tss": None} }, upsert = True)
			    #print "Found m: {0}".format(m)
			    mcollection.insert({"motif_id": wmid, "tf_name": tfName, "genomic_regions_gene_mapping": m["genomic_regions_gene_mapping"], "ct_info": m["ct_info"], "motif_tf_info": m["motif_tf_info"], "motif_genomic_regions_info": {"chr": wchrom, "start": wstart, "end": wend, "strand": wstrand}, "motif_score": wpval})  
			    mcollection.remove({"motif_id": row[0], "motif_score": None})
			    #print mcollection.count()
			    #ofile.writelines('\t'.join(modRow(wrow))+'\n')
		    else:
			if len(rowlist) > 1:
			    wrow = rowlist[minP(rowlist)]
			    ##best of the overlapping motifs, current line non-overlapping
			elif len(rowlist) == 1:
			    wrow = prow
			    ##previous and current line non-overlapping 
			rowlist = []##reset container
			rowlist.append(row)
			wmid, wchrom, wstart, wend, wstrand, wpval = wrow[0], wrow[1], int(wrow[2]), int(wrow[3]), wrow[4], float(wrow[6])
			m = mcollection.find_one({"motif_id": wmid, "tf_name": tfName})
			mcollection.insert({"motif_id": wmid, "tf_name": tfName, "genomic_regions_gene_mapping": m["genomic_regions_gene_mapping"], "ct_info": m["ct_info"], "motif_tf_info": m["motif_tf_info"], "motif_genomic_regions_info": {"chr": wchrom, "start": wstart, "end": wend, "strand": wstrand}, "motif_score": wpval})
			#print mcollection.count()
			#mcollection.update({"motif_id": wmid, "tf_name": tfName}, {"$set": {"motif_genomic_regions_info.chr": wchrom, "motif_genomic_regions_info.start": start, "motif_genomic_regions_info.end": end, "motif_genomic_regions_info.strand": strand, "motif_score": wpval}})
			#ofile.writelines('\t'.join(modRow(wrow))+'\n')
		    #mcollection.remove({"motif_id": row[0], "motif_score": None})
		prow = row
		psize = size
		if len(rowlist) == 0:
		    rowlist.append(prow)##push first row in, during looping this should never be empty
		pchrom, pstart, pend, pstrand, ppval = chrom, start, end, strand, pval
	        #remove original placeholder entry
    cursor = mcollection.find()
    print 'updated count ', mcollection.count()
    #for item in cursor:
	#mcollection.insert({motif_id: item['motif_id']})
	#print item
    #    print item['tf_name']
	#c = cursor.next()
	#print cursor.next()	
	#collection.find_one({"motif_id":}) 

##update motif_genomic_regions_info
    #for 
	#db.motif_instance_hughes_test.insert({"motif_tf_info": db.motif_instance_hughes_test.findOne({'tf_name':'Irf4'}).motif_tf_info})

if __name__=='__main__':
    sys.exit(main(sys.argv))
