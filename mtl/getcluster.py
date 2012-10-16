import sys
import os
import csv
#import operator

def main(argv):
    if len(argv) < 2:
        sys.stderr.write("Usage: %s MTLS_File Cluster_Assignments\n" % argv[0])
        return 1
    if not os.path.isfile(argv[1]):
        sys.stderr.write('Error: MTLS_file %r was not found!\n' % argv[1])
        return 1
    if not os.path.isfile(argv[2]):
        sys.stderr.write('Error: Cluster_Assignments %r was not found!\n' % argv[2])
        return 1
##    if not os.path.isfile(argv[3]):
##        sys.stderr.write('Error: nm#_conversion_file %r was not found!\n' % argv[3])
##        return 1
    
    infile = sys.argv[1]
    infile2= sys.argv[2]

    (path,fname) = os.path.split(infile)
    (shortname, extension) = os.path.splitext(fname)
    ifile = open(infile,'rU')
    ifile2= open(infile2,'rU')
    reader = csv.reader(ifile, delimiter = '\t')
    reader2=csv.reader(ifile2, delimiter = '\t')
    ##writer = csv.writer(ofile, delimiter = '\t')

    ##(path1,fname1) = os.path.split(infile1)
    ##(path2,fname2) = os.path.split(infile2)

##ifile = open('/home/xc406/data/fimo052912.txt','rt')
##reader = csv.reader(ifile, delimiter = '\t')

    ofile = open('/home/xc406/data/mtl/' + 'cluster1', 'wt')
    writer = csv.writer(ofile, delimiter = '\t')
	
#sortedlist = sorted(reader, key=operator.itemgetter(3), reverse=False)	    

    ##write row_number/target(mtl.id) to a list
    cluster1 = []
    ifile.seek(0)
    for row in reader2:
	if row[1] == 'cluster_1':
	    cluster1.append(row[0])	
    
    print cluster1

    for row in reader:
	#if not row == pline:
	if row[0] in cluster1:
	    newrow = ['','','']
	    newrow[0] = row[1]
	    newrow[1] = str(int(row[4])-100)
	    newrow[2] = str(int(row[5])+100)
            #newrow[3] = row[3]
	    #newrow[4] = row[0]  
            writer.writerows([newrow])


    ofile.close()

if __name__=='__main__':
    sys.exit(main(sys.argv))



