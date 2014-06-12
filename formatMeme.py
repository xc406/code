import sys
import os
import csv

def main(argv):
    if len(argv) < 3:
        sys.stderr.write("Usage: %s meme_input_file TF_Info_file \n" % argv[0])
        return 1
    if not os.path.isfile(argv[1]):
        sys.stderr.write('Error: meme_input_file %r was not found!\n' % argv[1])
        return 1
    if not os.path.isfile(argv[2]):
        sys.stderr.write('Error: TF_Info_file %r was not found!\n' % argv[2])
        return 1
    #if not os.path.isfile(argv[3]):
        #sys.stderr.write('Error: nm#_conversion_file %r was not found!\n' % argv[3])
        #return 1

    ##write meme/fimo output to fimoutall by column duplicating 1st and 2nd columns 
    infile = sys.argv[1]

    (path,fname) = os.path.split(infile)

    ifile = open(infile,'rt')
    ofile = open(os.path.join(path, "mm9motifs90cutmaptf.meme"),'w')
    reader = csv.reader(ifile, delimiter = '\t')
    writer = csv.writer(ofile, delimiter = '\t')

    #writer.writerows([[row[i] for i in [0,0,1,1,2,3,-4,-3,-2,-1]] for row in reader])

    ##substitute motif_id and nm# with tf and gene_name respectively
    infile1 = sys.argv[2]

    ##(path1,fname1) = os.path.split(infile1)
    ##(path2,fname2) = os.path.split(infile2)

    ifile1 = open(infile1,'rt')

    reader1 = csv.reader(ifile1, delimiter = '\t')

    ## generate a list of tuples (motif_id,tf)
    mylist1 = list([(row[3],row[6]) for row in reader1]) 
    #print mylist1

    ##store motif_id, tf in mydict1
    mydict1 = {}
    for k,v in mylist1:
        if not k in mydict1:
                mydict1[k] = v
        else:
                mydict1[k] += ','
                mydict1[k] += v

    keylist = mydict1.keys()
    keylist.sort()
    d1 = mydict1.items()    
    
    print mydict1
    ##substitute motif_id with tf
    ifile.seek(0)
    count = 0
    for row in reader:
	print row
	if len(row) >0:
	    if 'MOTIF' in row[0]:
	        rlist = row[0].split(' ')
	        #print 'before', rlist
                if rlist[0] == 'MOTIF':
	            for k,v in d1:
                        if k == rlist[1]:
			    count += 1
                            rlist[2] = v
		            #print 'after', rlist
		rstr = (' ').join(rlist)
                ofile.write(rstr+'\n')
	    else:
		rstr = (' ').join(row)
	        ofile.write(rstr+'\n')
	else:
	    ofile.write('\n')
    print count
    ofile.close()

if __name__=='__main__':
    sys.exit(main(sys.argv))
   
