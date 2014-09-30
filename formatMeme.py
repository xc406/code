"""format meme files with [MOTIF motif_id tf_name]""" 
import sys, os, csv, glob

def main(argv):
    if len(argv) < 3:
        sys.stderr.write("Usage: %s meme_input_file TF_Info_file \n" % argv[0])
        return 1
    if not os.path.exists(argv[1]):
        sys.stderr.write('Error: meme_input_file_path %r was not found!\n' % argv[1])
        return 1
    if not os.path.isfile(argv[2]):
        sys.stderr.write('Error: TF_Info_file %r was not found!\n' % argv[2])
        return 1
    #if not os.path.isfile(argv[3]):
        #sys.stderr.write('Error: nm#_conversion_file %r was not found!\n' % argv[3])
        #return 1

    ##write meme/fimo output to fimoutall by column duplicating 1st and 2nd columns 
    path = sys.argv[1]

    #(path,fname) = os.path.split(infile)

    #ifile = open(infile,'rt')
    #ofile = open(os.path.join(path, "mm9motifs90cutmaptf.meme"),'w')
    #reader = csv.reader(ifile, delimiter = '\t')
    #writer = csv.writer(ofile, delimiter = '\t')

    #writer.writerows([[row[i] for i in [0,0,1,1,2,3,-4,-3,-2,-1]] for row in reader])
	
    ##substitute motif_id and nm# with tf and gene_name respectively
    tfInfoFile = sys.argv[2]

    ##(path1,fname1) = os.path.split(infile1)
    ##(path2,fname2) = os.path.split(infile2)

    tfInfoReader = csv.reader(open(tfInfoFile,'rt'), delimiter = '\t')

    ## generate a list of tuples (motif_id,tf)
    mylist1 = list([(row[3],row[6]) for row in tfInfoReader]) 
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

    for infile in glob.glob(os.path.join(path,"*.meme")):
	fname = os.path.split(infile)[1]
	reader = csv.reader(open(infile,'rt'), delimiter = '\t')
	ofile = open(os.path.join(path+"/mapped_meme",fname),'wt')
	for row in reader:
	    if len(row) > 0 and 'MOTIF' in row[0]:
		rlist = row[0].split(' ')
		if '_1.00' in rlist[1]:
		    rlist[2] = mydict1[rlist[1]]
		elif '>' in rlist[1]:
		    info = row[0].split('>')[-1]
		    rlist = [rlist[0],rlist[3],'>'+info]
		else:
		    print "Warning:", row
	        rstr = (' ').join(rlist)
	        ofile.write(rstr+'\n')
	    else:
		rstr = (' ').join(row)
		ofile.write(rstr+'\n')	
	ofile.close()
    ##substitute motif_id with tf
    #ifile.seek(0)
    #count = 0
    #for row in reader:
#	print row
#	if len(row) >0:
#	    if 'MOTIF' in row[0]:
#	        rlist = row[0].split(' ')
#	        #print 'before', rlist
#                if rlist[0] == 'MOTIF':
#	            for k,v in d1:
#                        if k == rlist[1]:
#			    count += 1
#                            rlist[2] = v
#		            #print 'after', rlist
#		rstr = (' ').join(rlist)
#                ofile.write(rstr+'\n')
#	    else:
#		rstr = (' ').join(row)
#	        ofile.write(rstr+'\n')
#	else:
#	    ofile.write('\n')
#    print count
#    ofile.close()

if __name__=='__main__':
    sys.exit(main(sys.argv))
   
