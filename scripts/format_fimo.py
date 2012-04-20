import sys
import os
import csv

def main(argv):
    if len(argv) < 4:
        sys.stderr.write("Usage: %s fimo_output_file TF_Info_file nm#_conversion_file\n" % argv[0])
        return 1
    if not os.path.isfile(argv[1]):
        sys.stderr.write('Error: fimo_output_file %r was not found!\n' % argv[1])
        return 1
    if not os.path.isfile(argv[2]):
        sys.stderr.write('Error: TF_Info_file %r was not found!\n' % argv[2])
        return 1
    if not os.path.isfile(argv[3]):
        sys.stderr.write('Error: nm#_conversion_file %r was not found!\n' % argv[3])
        return 1

    ##write meme/fimo output to fimoutall by column duplicating 1st and 2nd columns 
    infile = sys.argv[1]

    (path,fname) = os.path.split(infile)

    ifile = open(infile,'rt')
    ofile = open(os.path.join(path, "fimoutall"),'w')
    reader = csv.reader(ifile, delimiter = '\t')
    writer = csv.writer(ofile, delimiter = '\t')

    writer.writerows([[row[i] for i in [0,0,1,1,2,3,-4,-3,-2,-1]] for row in reader])

    ifile.close()
    ofile.close()

    ##substitute motif_id and nm# with tf and gene_name respectively
    infile1 = sys.argv[2]
    infile2 = sys.argv[3]

    ##(path1,fname1) = os.path.split(infile1)
    ##(path2,fname2) = os.path.split(infile2)

    ifile1 = open(infile1,'rt')
    ifile2 = open(infile2,'rt')

    ofileo = open(os.path.join(path, "hg19"), 'w')
    writero = csv.writer(ofileo, delimiter = '\t')

    ofile = open(os.path.join(path, "fimo_formatted_all"),'r')
    readero = csv.reader(ofile, delimiter = '\t')
    reader1 = csv.reader(ifile1, delimiter = '\t')
    reader2 = csv.reader(ifile2, delimiter = '\t')

    a=0
    b=0
    c=0
    i=0
    j=0
    x=0
    y=0
    z=0

    ## generate a list of tuples (motif_id,tf)
    mylist1 = list([(row[3],row[6]) for row in reader1]) 
    ##print mylist1

    ## generate a list of tuples (nm#,gene_name)
    mylist2 = []
    mylist2 = list([(row[1],row[-4]) for row in reader2]) 
    ##print mylist2
    ##print len(mylist2)

    ##store motif_id, tf in mydict1
    mydict1 = {}
    for k,v in mylist1:
        if not k in mydict1:
                mydict1[k] = v
        else:
                mydict1[k] += ', '
                mydict1[k] += v

    keylist = mydict1.keys()
    keylist.sort()
    d1 = mydict1.items()    

    ##store nm#, gene_names in mydict2  
    mydict2 = {}
    for k,v in mylist2:
        if not k in mydict2:
                mydict2[k] = v
    ##	else:
    ##		mydict2[k] += ', '
    ##		mydict2[k] += v##duplicates

    d2 = mydict2.items()

    ##substitute motif_id with tf
    ##substitute nm# with gene_names
    ofile.seek(0)
    for row in readero:
        for k,v in d1:
            if k == row[0]:
                row[1] = v
                a+=1
            b+=1

        for k,v in d2:
            if row[2] == k:
                row[3] = v
                x+=1
            y+=1
        writero.writerows([row])
        z+=1

    print a, b, x, y, z

    ofileo.close()

if __name__=='__main__':
    sys.exit(main(sys.argv))
   
