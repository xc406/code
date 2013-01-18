import sys
import os
import csv

def main(argv):
    if len(argv) < 3:
        sys.stderr.write("Usage: %s sam_diff_expr_file mgd_file\n" % argv[0])
        return 1
    if not os.path.isfile(argv[1]):
        sys.stderr.write('Error: sam_diff_expr_file %r was not found!\n' % argv[1])
        return 1
    if not os.path.isfile(argv[2]):
        sys.stderr.write('Error: mgd_file %r was not found!\n' % argv[2])
        return 1
    #if not os.path.isfile(argv[3]):
        #sys.stderr.write('Error: nm#_conversion_file %r was not found!\n' % argv[3])
        #return 1

    ##write meme/fimo output to fimoutall by column duplicating 1st and 2nd columns 
    infile = sys.argv[1]

    (path,fname) = os.path.split(infile)
    (name,ext) = os.path.splitext(fname)

    ifile = open(infile,'U')
    ofile = open(os.path.join(path, name+"mm.txt"),'w')
    reader = csv.reader(ifile, delimiter = '\t')
    writer = csv.writer(ofile, delimiter = '\t')

    ##substitute hg with mm gene_names
    infile1 = sys.argv[2] ## conversion file from MGD
    
    ##(path1,fname1) = os.path.split(infile1)
    ##(path2,fname2) = os.path.split(infile2)

    ifile1 = open(infile1,'rt')
    #ifile2 = open(infile2,'rt')

    #ofile = open(os.path.join(path, "fimoutall"),'r')
    #readero = csv.reader(ofile, delimiter = '\t')
    reader1 = csv.reader(ifile1, delimiter = '\t')
    #reader2 = csv.reader(ifile2, delimiter = '\t')

    ## generate a list of tuples (hg,mm)
    #mylist1 = list([(row[3],row[0]) for row in reader1]) 
    ##(mm,hg)
    mylist1 = list([(row[0],row[3]) for row in reader1])
    
    ##store {hg_name, mm_name} in mydict1
    mydict1 = {}
    for k,v in mylist1:
        if not k.upper() in mydict1:
                mydict1[k.upper()] = v
    keylist = mydict1.keys()
    keylist.sort()
    d1 = mydict1.items()    
    #print len(keylist),len(mylist1)
    #print mylist1

    for row in reader:
	#print row[0]
	if row[0].upper() in keylist:
		#print row[0]
		#row[0] = mydict1[row[0].upper()]
		row[0] = row[0].upper() 
	#if ('_0.62' in row[0]) and (row[0] in keylist):
        #    if ',' in mydict1[row[0]]:
	#	for i in range(len(mydict1[row[0]].split(', '))):
	#	    if not mydict1[row[0]].split(', ')[i] in mydict2:
	#	        mydict2[mydict1[row[0]].split(', ')[i]] = float(row[1])
	#	    else:
	#		mydict2[mydict1[row[0]].split(', ')[i]] = min(float(row[1]), mydict2[mydict1[row[0]].split(', ')[i]])
	#    else:
	#	if not mydict1[row[0]] in mydict2:
	#	    mydict2[mydict1[row[0]]] = float(row[1])
	#	else:
	#	    mydict2[mydict1[row[0]]] = min(float(row[1]), mydict2[mydict1[row[0]]])
		writer.writerows([row])
    ofile.close()

if __name__=='__main__':
    sys.exit(main(sys.argv))
   
