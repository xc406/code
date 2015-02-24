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
                item.get("cons").get("rSnp",'NA'),
                item.get("cons").get("rIndel",'NA'),
                item.get("cons").get("rphastCons100way",'NA'),
                item.get("cons").get("phyloP100way",'NA'),
                item.get("cons").get("rphastCons46way_primates",'NA'),
                item.get("cons").get("phyloP46way_primate",'NA'),
		item.get("cons").get("rphastCons46way_placental",'NA'),
		item.get("cons").get("phyloP46way_placental",'NA'),
		item.get("cons").get("rphastCons46way",'NA'),
		item.get("cons").get("phyloP46way",'NA'),
		item.get("cons").get("piSnp",'NA'),
                item.get("cons").get("piIndel",'NA'),
                item.get("cons").get("phastCons100way",'NA'),
                item.get("cons").get("phastCons46way_primates",'NA'),
                item.get("cons").get("phastCons46way_placental",'NA'),
                item.get("cons").get("phastCons46way",'NA'),
                item["gc"][0],
                item["gc"][1],
                #item["map"]["wgEncodeDukeMapabilityUniqueness35bpchr17"],
                #item["map"]["wgEncodeDukeMapabilityUniqueness20bpchr17"],
                #item["map"]["wgEncodeCrgMapabilityAlign36merchr17"],
                item["map"]["wgEncodeCrgMapabilityAlign24mer"],
                item["dgf"].values()[1],#["A549Aln"],#.values()[1],#["K562"],
		item.get("dgf").get("fos","NA"),#"A549fos","NA"),
		item.get("dgf").get("fosStrand",'NA'),#"A549fosStrand",'NA'),
		#item["dnase"]["Gm12878AlnRep1"]+item["dnase"]["Gm12878AlnRep1"],
		#item.get("dnase").get("fos","NA"),
		#item.get("dnase").get("fosStrand","NA"),
		#item["dnase"]["Gm12878AlnRep1"],
		#item["dnase"]["Gm12878AlnRep2"],
                "".join(item.get("map").get("wgEncodeDukeMapabilityRegionsExcludable",'NA')),
                "".join(item.get("map").get("wgEncodeDacMapabilityConsensusExcludable",'NA')),
		#item["motif_id"].split(".")[0],
		#item.get("exon","NA"),
		"".join(item.get("exon","NA")),
		len(item.get("target_gene").get("10kb",[])),
		"NA","NA","NA","NA","NA","NA","NA","NA","NA","NA",
		"NA","NA","NA","NA","NA","NA","NA","NA","NA","NA",
		"NA","NA","NA","NA","NA","NA","NA","NA","NA"#,"NA",
		#"NA","NA","NA","NA","NA","NA","NA","NA","NA","NA",
		#"NA","NA"
		##29 for K549; 35(29+6treated) for A549; 42 for Gm12878
                ]
	else:
		row=[item["genomic_region"]["chr"],
		item["genomic_region"]["start"],
		item["genomic_region"]["end"],
		#item["genomic_region"]["strand"],
		item["motif_score"],
		item["target_gene"]["dist_tss"],
		item.get("cons").get("rSnp",'NA'),
		item.get("cons").get("rIndel",'NA'),
		item.get("cons").get("rphastCons100way",'NA'),
		item.get("cons").get("phyloP100way",'NA'),
		item.get("cons").get("rphastCons46way_primates",'NA'),
		item.get("cons").get("phyloP46way_primate",'NA'),
		item.get("cons").get("rphastCons46way_placental",'NA'),
		item.get("cons").get("phyloP46way_placental",'NA'),
		item.get("cons").get("rphastCons46way",'NA'),
		item.get("cons").get("phyloP46way",'NA'),
		item.get("cons").get("piSnp",'NA'),
                item.get("cons").get("piIndel",'NA'),
                item.get("cons").get("phastCons100way",'NA'),
                item.get("cons").get("phastCons46way_primates",'NA'),
                item.get("cons").get("phastCons46way_placental",'NA'),
                item.get("cons").get("phastCons46way",'NA'),
		item["gc"][0],
		item["gc"][1],
		#item["map"]["wgEncodeDukeMapabilityUniqueness35bpchr17"],
		#item["map"]["wgEncodeDukeMapabilityUniqueness20bpchr17"],
		#item["map"]["wgEncodeCrgMapabilityAlign36merchr17"],
		item["map"]["wgEncodeCrgMapabilityAlign24mer"],
		item["dgf"].values()[1],#["K562"],#["A549Aln"],
		item.get("dgf").get("fos","NA"),
		item.get("dgf").get("fosStrand",'NA'),#"A549fosStrand",'NA'),
		#item["dnase"]["Gm12878AlnRep1"]+item["dnase"]["Gm12878AlnRep1"],
		#item.get("dnase").get("fos","NA"),
		#item.get("dnase").get("fosStrand",'NA'),
		#item["dnase"]["Gm12878AlnRep1"],
		#item["dnase"]["Gm12878AlnRep2"],
		"".join(item.get("map").get("wgEncodeDukeMapabilityRegionsExcludable",'NA')),
		"".join(item.get("map").get("wgEncodeDacMapabilityConsensusExcludable",'NA')),
		#item["motif_id"].split('.')[0],
		#item.get("exon","NA"),
		"".join(item.get("exon","NA")),
		len(item.get("target_gene").get("10kb",[])),
		#item.get("chip").get("wgEncodeAwgTfbsHaibGm12878Yy1sc281Pcr1xUniPkNarrowPeakFormat","NA"),
		#item.get("chip").get("wgEncodeAwgTfbsSydhGm12878Yy1UniPkNarrowPeakFormat","NA"),#caution
		#item.get("chip").get("wgEncodeHaibTfbsGm12878Yy1sc281Pcr1xPkRep1Format","NA"),
		#item.get("chip").get("wgEncodeHaibTfbsGm12878Yy1sc281Pcr1xPkRep2Format","NA"),
		#item.get("chip").get("wgEncodeSydhTfbsGm12878Yy1StdPkNarrowPeakFormat","NA"),
		#item.get("chip").get("wgEncodeAwgTfbsSydhGm12878Znf143166181apUniPkNarrowPeakFormat","NA"),
		#item.get("chip").get("wgEncodeSydhTfbsGm12878Znf143166181apStdPkNarrowPeakFormat","NA"),
		#item.get("chip").get("wgEncodeSydhTfbsGm12878Znf384hpa004051IggmusPkNarrowPeakFormat","NA"),
		#item.get("chip").get("wgEncodeAwgTfbsSydhGm12878NfyaIggmusUniPkNarrowPeakFormat","NA"),
		#item.get("chip").get("wgEncodeSydhTfbsGm12878NfyaIggmusPkNarrowPeakFormat","NA"),
		#item.get("chip").get("wgEncodeAwgTfbsSydhGm12878E2f4IggmusUniPkNarrowPeakFormat","NA"),
		#item.get("chip").get("wgEncodeSydhTfbsGm12878E2f4IggmusPkNarrowPeakFormat","NA"),
		#item.get("chip").get("wgEncodeAwgTfbsSydhGm12878JundUniPkNarrowPeakFormat","NA"),
		#item.get("chip").get("wgEncodeSydhTfbsGm12878JundIggrabPkNarrowPeakFormat","NA"),
		#item.get("chip").get("wgEncodeSydhTfbsGm12878JundStdPkNarrowPeakFormat","NA"),
		#item.get("chip").get("wgEncodeBroadHistoneGm12878CtcfStdPkFormat","NA"),
		#item.get("chip").get("wgEncodeUwTfbsGm12878CtcfStdPkRep1NarrowPeakFormat","NA"),
		#item.get("chip").get("wgEncodeUwTfbsGm12878CtcfStdPkRep2NarrowPeakFormat","NA"),
		#item.get("chip").get("wgEncodeAwgTfbsBroadGm12878CtcfUniPkNarrowPeakFormat","NA"),
		#item.get("chip").get("wgEncodeAwgTfbsSydhGm12878Ctcfsc15914c20UniPkNarrowPeakFormat","NA"),
		#item.get("chip").get("wgEncodeAwgTfbsUtaGm12878CtcfUniPkNarrowPeakFormat","NA"),
		#item.get("chip").get("wgEncodeAwgTfbsUwGm12878CtcfUniPkNarrowPeakFormat","NA"),
		#item.get("chip").get("wgEncodeSydhTfbsGm12878Ctcfsc15914c20StdPkNarrowPeakFormat","NA"),
		#item.get("chip").get("wgEncodeUwTfbsGm12878CtcfStdHotspotsRep1Format","NA"),
		#item.get("chip").get("wgEncodeUwTfbsGm12878CtcfStdHotspotsRep2Format","NA"),
		#item.get("chip").get("wgEncodeSydhTfbsGm12878MafkIggmusPkNarrowPeakFormat","NA"),
		#item.get("chip").get("wgEncodeAwgTfbsSydhGm12878MaxIggmusUniPkNarrowPeakFormat","NA"),
		#item.get("chip").get("wgEncodeSydhTfbsGm12878MaxIggmusPkNarrowPeakFormat","NA"),
		#item.get("chip").get("wgEncodeSydhTfbsGm12878MaxStdPkNarrowPeakFormat","NA"),
		#item.get("chip").get("wgEncodeAwgTfbsHaibGm12878NrsfPcr1xUniPkNarrowPeakFormat","NA"),
		#item.get("chip").get("wgEncodeHaibTfbsGm12878NrsfPcr1xPkRep1Format","NA"),
		#item.get("chip").get("wgEncodeHaibTfbsGm12878NrsfPcr1xPkRep2Format","NA"),
		#item.get("chip").get("wgEncodeHaibTfbsGm12878NrsfPcr2xPkRep1Format","NA"),
		#item.get("chip").get("wgEncodeHaibTfbsGm12878NrsfPcr2xPkRep2Format","NA"),
		#item.get("chip").get("wgEncodeAwgTfbsSydhGm12878Sin3anb6001263IggmusUniPkNarrowPeakFormat","NA"),
		#item.get("chip").get("wgEncodeSydhTfbsGm12878Sin3anb6001263IggmusPkNarrowPeakFormat","NA"),
		#item.get("chip").get("wgEncodeAwgTfbsHaibGm12878Sp1Pcr1xUniPkNarrowPeakFormat","NA"),
		#item.get("chip").get("wgEncodeHaibTfbsGm12878Sp1Pcr1xPkRep1Format","NA"),
		#item.get("chip").get("wgEncodeHaibTfbsGm12878Sp1Pcr1xPkRep2Format","NA"),
		#item.get("chip").get("wgEncodeAwgTfbsHaibGm12878Usf1Pcr2xUniPkNarrowPeakFormat","NA"),
		#item.get("chip").get("wgEncodeHaibTfbsGm12878Usf1Pcr2xPkRep1Format","NA"),
		#item.get("chip").get("wgEncodeHaibTfbsGm12878Usf1Pcr2xPkRep2Format","NA")

		#item.get("chip").get("wgEncodeAwgTfbsHaibA549Ctcfsc5916Pcr1xEtoh02UniPk","NA"),#"wgEncodeAwgTfbsUwK562CtcfUniPk",'NA'),
		#item.get("chip").get("wgEncodeAwgTfbsUtaA549CtcfUniPk","NA"),#"wgEncodeAwgTfbsSydhK562CtcfbIggrabUniPk",'NA'),
		#item.get("chip").get("wgEncodeAwgTfbsUwA549CtcfUniPk","NA"),#wgEncodeAwgTfbsBroadK562CtcfUniPk",'NA'),
		#item.get("chip").get("wgEncodeHaibTfbsA549Ctcfsc5916Pcr1xEtoh02PkRep1","NA"),#"wgEncodeAwgTfbsHaibK562CtcfcPcr1xUniPk",'NA'),
		#item.get("chip").get("wgEncodeHaibTfbsA549Ctcfsc5916Pcr1xEtoh02PkRep2","NA"),#"wgEncodeAwgTfbsUtaK562CtcfUniPk",'NA'),
		#item.get("chip").get("wgEncodeSydhTfbsA549CtcfbIggrabPk","NA"),#"wgEncodeSydhTfbsK562CtcfbIggrabPk",'NA'),
		#item.get("chip").get("wgEncodeUwTfbsA549CtcfStdPkRep1","NA"),#"wgEncodeAwgTfbsSydhK562E2f4UcdUniPk","NA"),
		#item.get("chip").get("wgEncodeUwTfbsA549CtcfStdPkRep2","NA"),#"wgEncodeSydhTfbsK562E2f4UcdPk","NA"),
		#item.get("chip").get("wgEncodeHaibTfbsA549JundV0416102Etoh02PkRep1","NA"),#"wgEncodeAwgTfbsSydhK562JundIggrabUniPk","NA"),
		#item.get("chip").get("wgEncodeHaibTfbsA549JundV0416102Etoh02PkRep2","NA"),#"wgEncodeSydhTfbsK562JundIggrabPk","NA"),
		#item.get("chip").get("wgEncodeAwgTfbsSydhA549MaxIggrabUniPk","NA"),#"wgEncodeAwgTfbsSydhK562Mafkab50322IggrabUniPk","NA"),
		#item.get("chip").get("wgEncodeHaibTfbsA549MaxV0422111PkRep1","NA"),#"wgEncodeSydhTfbsK562Mafkab50322IggrabPk","NA"),
		#item.get("chip").get("wgEncodeHaibTfbsA549MaxV0422111PkRep2","NA"),#"wgEncodeAwgTfbsHaibK562MaxV0416102UniPk","NA"),
		#item.get("chip").get("wgEncodeSydhTfbsA549MaxIggrabPk","NA"),#"wgEncodeSydhTfbsK562MaxIggrabPk","NA"),
		#item.get("chip").get("wgEncodeAwgTfbsHaibA549NrsfV0422111Etoh02UniPk","NA"),#"wgEncodeAwgTfbsSydhK562MaxIggrabUniPk","NA"),
		#item.get("chip").get("wgEncodeHaibTfbsA549NrsfV0422111Etoh02PkRep1","NA"),#"wgEncodeSydhTfbsK562MaxStdPk","NA"),
		#item.get("chip").get("wgEncodeHaibTfbsA549NrsfV0422111Etoh02PkRep2","NA"),#"wgEncodeAwgTfbsSydhK562NfyaUniPk","NA"),
		#item.get("chip").get("wgEncodeAwgTfbsHaibA549Sin3ak20V0422111Etoh02UniPk","NA"),#"wgEncodeSydhTfbsK562NfyaStdPk","NA"),
		#item.get("chip").get("wgEncodeHaibTfbsA549Sp1V0422111Etoh02PkRep1","NA"),#"wgEncodeAwgTfbsHaibK562NrsfV0416102UniPk","NA"),
		#item.get("chip").get("wgEncodeHaibTfbsA549Sp1V0422111Etoh02PkRep2","NA"),#"wgEncodeAwgTfbsHaibK562Sin3ak20V0416101UniPk","NA"),
		#item.get("chip").get("wgEncodeAwgTfbsHaibA549Usf1Pcr1xEtoh02UniPk","NA"),#"wgEncodeAwgTfbsHaibK562Sp1Pcr1xUniPk","NA"),
		#item.get("chip").get("wgEncodeAwgTfbsHaibA549Usf1V0422111Etoh02UniPk","NA"),#"wgEncodeAwgTfbsHaibK562Usf1V0416101UniPk","NA"),
		#item.get("chip").get("wgEncodeHaibTfbsA549Usf1Pcr1xEtoh02PkRep1","NA"),#"wgEncodeAwgTfbsHaibK562Yy1V0416101UniPk","NA"),
		#item.get("chip").get("wgEncodeHaibTfbsA549Usf1Pcr1xEtoh02PkRep2","NA"),#"wgEncodeAwgTfbsSydhK562Yy1UcdUniPk","NA"),
		#item.get("chip").get("wgEncodeHaibTfbsA549Usf1V0422111Etoh02PkRep1","NA"),#"wgEncodeAwgTfbsHaibK562Yy1V0416102UniPk","NA"),
		#item.get("chip").get("wgEncodeHaibTfbsA549Usf1V0422111Etoh02PkRep2","NA"),#"wgEncodeSydhTfbsK562Yy1UcdPk","NA"),
		#item.get("chip").get("wgEncodeAwgTfbsHaibA549Yy1cV0422111Etoh02UniPk","NA"),#"wgEncodeAwgTfbsSydhK562Znf143IggrabUniPk","NA"),
		#item.get("chip").get("wgEncodeHaibTfbsA549Yy1cV0422111Etoh02PkRep1","NA"),#"wgEncodeSydhTfbsK562Znf143IggrabPk","NA"),		
		#item.get("chip").get("wgEncodeHaibTfbsA549Yy1cV0422111Etoh02PkRep2","NA")		

		item.get("chip").get("wgEncodeAwgTfbsUwK562CtcfUniPk",'NA'),
		item.get("chip").get("wgEncodeAwgTfbsSydhK562CtcfbIggrabUniPk",'NA'),
		item.get("chip").get("wgEncodeAwgTfbsBroadK562CtcfUniPk",'NA'),
		item.get("chip").get("wgEncodeAwgTfbsHaibK562CtcfcPcr1xUniPk",'NA'),
		item.get("chip").get("wgEncodeAwgTfbsUtaK562CtcfUniPk",'NA'),
		item.get("chip").get("wgEncodeSydhTfbsK562CtcfbIggrabPk",'NA'),
		item.get("chip").get("wgEncodeAwgTfbsSydhK562E2f4UcdUniPk","NA"),
		item.get("chip").get("wgEncodeSydhTfbsK562E2f4UcdPk","NA"),
		item.get("chip").get("wgEncodeAwgTfbsSydhK562JundIggrabUniPk","NA"),
		item.get("chip").get("wgEncodeSydhTfbsK562JundIggrabPk","NA"),
		item.get("chip").get("wgEncodeAwgTfbsSydhK562Mafkab50322IggrabUniPk","NA"),
		item.get("chip").get("wgEncodeSydhTfbsK562Mafkab50322IggrabPk","NA"),
		item.get("chip").get("wgEncodeAwgTfbsHaibK562MaxV0416102UniPk","NA"),
		item.get("chip").get("wgEncodeSydhTfbsK562MaxIggrabPk","NA"),
		item.get("chip").get("wgEncodeAwgTfbsSydhK562MaxIggrabUniPk","NA"),
		item.get("chip").get("wgEncodeSydhTfbsK562MaxStdPk","NA"),
		item.get("chip").get("wgEncodeAwgTfbsSydhK562NfyaUniPk","NA"),
		item.get("chip").get("wgEncodeSydhTfbsK562NfyaStdPk","NA"),
		item.get("chip").get("wgEncodeAwgTfbsHaibK562NrsfV0416102UniPk","NA"),
		item.get("chip").get("wgEncodeAwgTfbsHaibK562Sin3ak20V0416101UniPk","NA"),
		item.get("chip").get("wgEncodeAwgTfbsHaibK562Sp1Pcr1xUniPk","NA"),
		item.get("chip").get("wgEncodeAwgTfbsHaibK562Usf1V0416101UniPk","NA"),
		item.get("chip").get("wgEncodeAwgTfbsHaibK562Yy1V0416101UniPk","NA"),
		item.get("chip").get("wgEncodeAwgTfbsSydhK562Yy1UcdUniPk","NA"),
		item.get("chip").get("wgEncodeAwgTfbsHaibK562Yy1V0416102UniPk","NA"),
		item.get("chip").get("wgEncodeSydhTfbsK562Yy1UcdPk","NA"),
		item.get("chip").get("wgEncodeAwgTfbsSydhK562Znf143IggrabUniPk","NA"),
		item.get("chip").get("wgEncodeSydhTfbsK562Znf143IggrabPk","NA"),
		item.get("chip").get("wgEncodeSydhTfbsK562Znf384hpa004051IggrabPk","NA")
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
        row = [item["genomic_region"]["chr"],
                item["genomic_region"]["start"]-1-window,##0-based correction
                item["genomic_region"]["end"]+window,
                item["dgf"].values()[0]]
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
