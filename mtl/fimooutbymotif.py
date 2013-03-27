import sys
import os
import csv

def main(argv):
    if len(argv) < 2:
        sys.stderr.write("Usage: %s fimo_output_file\n" % argv[0])
        return 1
    if not os.path.isfile(argv[1]):
        sys.stderr.write('Error: fimo_output_file %r was not found!\n' % argv[1])
        return 1
##    if not os.path.isfile(argv[2]):
##        sys.stderr.write('Error: TF_Info_file %r was not found!\n' % argv[2])
##        return 1
##    if not os.path.isfile(argv[3]):
##        sys.stderr.write('Error: nm#_conversion_file %r was not found!\n' % argv[3])
##        return 1
    
    infile = sys.argv[1]

    (path,fname) = os.path.split(infile)
    (n,ext) = os.path.splitext(fname)
    ifile = open(infile,'rt')
    reader = csv.reader(ifile, delimiter = '\t')
    
    ofile = open('/home/xc406/data/mtl/' + n +'_count', 'wt')
    writer = csv.writer(ofile, delimiter = '\t')

    ##(path1,fname1) = os.path.split(infile1)
    ##(path2,fname2) = os.path.split(infile2)

    ##ifiletest = open('/home/xc406/data/mtl/cluster1fimoout_100812.txt','rt')
    ##readertest = csv.reader(ifiletf, delimiter = '\t')

    motiflist = []
    targetlist = []
    tempmo = ''
    motifcount = {}
    for row in reader:
	if '_0.71' in row[0]:
	    if row[0] != tempmo:
	        #motifcount[tempmo] = len(motiflist)
	        tempmo = row[0]
	        motiflist = []
		targetlist = []
	    elif row[0] == tempmo:
	        motiflist.append(row[0])
		if not row[1] in targetlist: 
		    targetlist.append(row[1])
		    motifcount[tempmo] = (len(motiflist),len(targetlist))
	            #motiflist.append(row[0])

    print motifcount

    for k in motifcount:
	newrow = [k,motifcount[k][1],motifcount[k][0]]
	writer.writerows([newrow])
    	
    ofile.close()
	
if __name__=='__main__':
    sys.exit(main(sys.argv))




