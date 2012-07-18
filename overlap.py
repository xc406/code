import sys
import os
import csv
#import operator

def main(argv):
    if len(argv) < 2:
        sys.stderr.write("Usage: %s gff_file\n" % argv[0])
        return 1
    if not os.path.isfile(argv[1]):
        sys.stderr.write('Error: gff_file %r was not found!\n' % argv[1])
        return 1
##    if not os.path.isfile(argv[2]):
##        sys.stderr.write('Error: TF_Info_file %r was not found!\n' % argv[2])
##        return 1
##    if not os.path.isfile(argv[3]):
##        sys.stderr.write('Error: nm#_conversion_file %r was not found!\n' % argv[3])
##        return 1
    
    infile = sys.argv[1]

    (path,fname) = os.path.split(infile)

    ifile = open(infile,'rt')
    reader = csv.reader(ifile, delimiter = '\t')
    ##writer = csv.writer(ofile, delimiter = '\t')

    ##(path1,fname1) = os.path.split(infile1)
    ##(path2,fname2) = os.path.split(infile2)

    ifilegn = open('/home/xc406/data/mm9_refseq_May_1_2011.txt','rt')
    readergn = csv.reader(ifilegn, delimiter = '\t')

    mylist1 = list([(row[1],row[-4]) for row in readergn])

    gndict = {}
    for m,n in mylist1:
        if not m in gndict:
            gndict[m]= n
        #else:
            #gndict[m] += ', '
            #gndict[m] += n

    nmlist = gndict.keys()

##ifile = open('/home/xc406/data/fimo052912.txt','rt')
##reader = csv.reader(ifile, delimiter = '\t')

    ofile = open('/home/xc406/data/mm9gff_example_gname/' +  fname, 'w')
    writer = csv.writer(ofile, delimiter = '\t')
    	
#		sortedlist = sorted(reader, key=operator.itemgetter(3), reverse=False)	    

    ##substitue NM# with hugo_gene_name
    ifile.seek(0)
    pline = ''
    for row in reader:
        #row[3] = str(int(row[3])-100)
	#row[4] = str(int(row[4])+100)
	group = row[-1].split(';')
	gid = group[0].split(' ')
	for k in nmlist:
	    if gid[1] == k:
		gid[1] = gndict[k]
		group[0] = gid[0] + ' ' + gid[1]
		group = group[0] + '; ' + 'transcript_id ' + k + ';' + group[1]
  		row[-1] = group
		if not row == pline: 
		    writer.writerows([row])
                    pline = row
    ofile.close()


if __name__=='__main__':
    sys.exit(main(sys.argv))



