import os, sys, csv
from array import array
from copy import copy
from collections import defaultdict
import time
from pymongo import MongoClient

def getData(cursor, writer, window=0):
    for item in cursor:
	if not 'chip' in item:
		row=[item["genomic_region"]["chr"],
                item["genomic_region"]["start"],
                item["genomic_region"]["end"],
                #item["genomic_region"]["strand"],
                item["motif_score"],
                item["target_gene"]["dist_tss"],
                item.get("cons").get("piSnp",'NA'),
                item.get("cons").get("piIndel",'NA'),
                item.get("cons").get("phastCons100way",'NA'),
                item.get("cons").get("phyloP100way",'NA'),
                item.get("cons").get("phastCons46way_primates",'NA'),
                item.get("cons").get("phyloP46way_primate",'NA'),
		item.get("cons").get("phastCons46way_placental",'NA'),
		item.get("cons").get("phyloP46way_placental",'NA'),
		item.get("cons").get("phastCons46way",'NA'),
		item.get("cons").get("phyloP46way",'NA'),
                item["gc"][0],
                item["gc"][1],
                item["map"]["wgEncodeDukeMapabilityUniqueness35bpchr17"],
                item["map"]["wgEncodeDukeMapabilityUniqueness20bpchr17"],
                item["map"]["wgEncodeCrgMapabilityAlign36merchr17"],
                item["map"]["wgEncodeCrgMapabilityAlign24merchr17"],
                item["Dgf"]["K562"],
		item.get("Dgf").get("fos",'NA'),
                'NA',
                'NA',
            	'NA',
                'NA',
                'NA',
                'NA',
                "".join(item.get("map").get("wgEncodeDukeMapabilityRegionsExcludable",'NA')),
                "".join(item.get("map").get("wgEncodeDacMapabilityConsensusExcludable",'NA')),
		item["motif_id"],
		item.get("exon","NA"),
		"".join(item.get("exon1st","NA"))
                ]
	else:
		row=[item["genomic_region"]["chr"],
		item["genomic_region"]["start"],
		item["genomic_region"]["end"],
		#item["genomic_region"]["strand"],
		item["motif_score"],
		item["target_gene"]["dist_tss"],
		item.get("cons").get("piSnp",'NA'),
		item.get("cons").get("piIndel",'NA'),
		item.get("cons").get("phastCons100way",'NA'),
		item.get("cons").get("phyloP100way",'NA'),
		item.get("cons").get("phastCons46way_primates",'NA'),
		item.get("cons").get("phyloP46way_primate",'NA'),
		item.get("cons").get("phastCons46way_placental",'NA'),
		item.get("cons").get("phyloP46way_placental",'NA'),
		item.get("cons").get("phastCons46way",'NA'),
		item.get("cons").get("phyloP46way",'NA'),
		item["gc"][0],
		item["gc"][1],
		item["map"]["wgEncodeDukeMapabilityUniqueness35bpchr17"],
		item["map"]["wgEncodeDukeMapabilityUniqueness20bpchr17"],
		item["map"]["wgEncodeCrgMapabilityAlign36merchr17"],
		item["map"]["wgEncodeCrgMapabilityAlign24merchr17"],
		item["Dgf"]["K562"],
		item.get("Dgf").get("fos",'NA'),
		item.get("chip").get("wgEncodeAwgTfbsUwK562CtcfUniPkchr17",'NA'),
		item.get("chip").get("wgEncodeAwgTfbsSydhK562CtcfbIggrabUniPkchr17",'NA'),
		item.get("chip").get("wgEncodeAwgTfbsBroadK562CtcfUniPkchr17",'NA'),
		item.get("chip").get("wgEncodeAwgTfbsHaibK562CtcfcPcr1xUniPkchr17",'NA'),
		item.get("chip").get("wgEncodeAwgTfbsUtaK562CtcfUniPkchr17",'NA'),
		item.get("chip").get("wgEncodeSydhTfbsK562CtcfbIggrabPkchr17",'NA'),
		"".join(item.get("map").get("wgEncodeDukeMapabilityRegionsExcludable",'NA')),
		"".join(item.get("map").get("wgEncodeDacMapabilityConsensusExcludable",'NA')),
		item["motif_id"],
		item.get("exon","NA"),
		"".join(item.get("exon1st","NA"))
		]
	writer.writerows([row])

def makeGff(cursor, gffWriter, window):
    """write gff files with a mongodb query"""
    for item in cursor:
        #print item["genomic_regions_gene_mapping"]
        info = 'gene_10kb ' + '_'.join(item["genomic_regions_gene_mapping"]["genelist10kb"]) + \
		'; closest_gene ' + '_'.join(item["genomic_regions_gene_mapping"]["closest_gene"])
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
	#print row
        gffWriter.writerows([row])

def makeBed4(cursor, bedWriter, window):
    """write bed files with a mongodb query"""
    for item in cursor:
        row = [item["motif_genomic_regions_info"]["chr"],
                item["motif_genomic_regions_info"]["start"]-1-window,##0-based correction
                item["motif_genomic_regions_info"]["end"]+window,
                item["motif_ct_info"]["accessibility_score"]]
        bedWriter.writerows([row])
	
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
	
	server = 'localhost'
	port = 27017
	client = MongoClient(server, port)
	db = client["mm9"]
	global mcollection
	mcollection = db["motif_instance_hughes_test"]
	startTime = time.clock()
	test = mcollection.find_one({"tf_name": "Zscan4"})
	window = 0
	test['new_field'] = {'new_new_field':'haha'}
	mcollection.save(test)
	#gff = open('/home/xc406/data/mongodbtest/testout2.gff','wt')
	#gffWriter = csv.writer(gff, delimiter = '\t')
	#makeGff(cursor, gffWriter, window) 
	print 'time', time.clock()-startTime
	# example usage
	#time1 = time.time()
	#wig = open(sys.argv[2],'rt')
	#wigFile = csv.reader(wig, delimiter = '\t')
	#ctName = "mapability"
	#coordDict, valuesDict = getBedCoord(wigFile, ctName)
	#dataType = 'phyloP46wayPrimate'
	
	#stepDict, startDict, valuesDict = getFixStart(wigFile, dataType)#'phyloP46wayPrimate')
	#dataType = 'phyloP46wayPrimate'
	#chrom = 'chr1'
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
