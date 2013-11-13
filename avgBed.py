import sys
import os
import csv
import numpy as np
"""get the average size of peaks(bed intervals)"""
def main(argv):
    if len(argv) < 2:
        sys.stderr.write("Usage: %s bed_file\n" % argv[0])
        return 1
    if not os.path.isfile(argv[1]):
        sys.stderr.write('Error: bed_file %r was not found!\n' % argv[1])
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

    #ifile = open('/home/xc406/data/mm9upstream10kb.bed','rt')
    #reader = csv.reader(ifile, delimiter = '\t')

    #ofile = open('/home/xc406/data/mm9ud2.gff', 'w')
    #writer = csv.writer(ofile, delimiter = '\t')
    
    a = np.array([])
    reader.next()
    for row in reader:
	size = int(row[2])-int(row[1])
	a = np.append(a,[size])
	#row.append('.')
	#row.append('gene_id '+ str(row[3]))
	#row.append('.')
	#row.append('.')
	#row[1],row[2],row[3],row[4],row[5] = row[3].split('_')[3],'fire',row[1],row[2],'1'
	#row[1],row[2],row[3],row[4],row[5],row[6],row[7],row[8] = str(row[0]+'_'+row[1]+'_'+row[2]),'TSS',row[1],row[2],'1','+','.',row[6]
	
	#writer.writerows([row])
    print len(a), np.mean(a)
    #ofile.close()


if __name__=='__main__':
    sys.exit(main(sys.argv))




