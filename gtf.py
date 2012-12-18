##path = '/Users/xclair/Documents/google-python-exercises/hg19'

import sys
import os
import csv

##def main(argv):
##    if len(argv) < 4:
##        sys.stderr.write("Usage: %s fimo_output_file TF_Info_file nm#_conversion_file\n" % argv[0])
##        return 1
##    if not os.path.isfile(argv[1]):
##        sys.stderr.write('Error: fimo_output_file %r was not found!\n' % argv[1])
##        return 1
##    if not os.path.isfile(argv[2]):
##        sys.stderr.write('Error: TF_Info_file %r was not found!\n' % argv[2])
##        return 1
##    if not os.path.isfile(argv[3]):
##        sys.stderr.write('Error: nm#_conversion_file %r was not found!\n' % argv[3])
##        return 1
##
##    ##write meme/fimo output to fimoutall by column duplicating 1st and 2nd columns 
##    infile = sys.argv[1]
##
##    (path,fname) = os.path.split(infile)
##
##    ifile = open(infile,'rt')
##    ofile = open(os.path.join(path, "fimoutall"),'w')
##    reader = csv.reader(ifile, delimiter = '\t')
##    writer = csv.writer(ofile, delimiter = '\t')
##
##    writer.writerows([[row[i] for i in [1,0,0,1,2,1,3,-4,-3,-2,-1]] for row in reader])
##
##    ifile.close()
##    ofile.close()
##
##    ##substitute motif_id and nm# with tf and gene_name respectively
##    infile1 = sys.argv[2]
##    infile2 = sys.argv[3]
##
##    ##(path1,fname1) = os.path.split(infile1)
##    ##(path2,fname2) = os.path.split(infile2)
##
##    ifile1 = open(infile1,'rt')
##    ifile2 = open(infile2,'rt')
##
##    ofileo = open(os.path.join(path, "hg19"), 'w')
##    writero = csv.writer(ofileo, delimiter = '\t')
##
##    ofile = open(os.path.join(path, "fimoutall"),'r')
##    readero = csv.reader(ofile, delimiter = '\t')
##    reader1 = csv.reader(ifile1, delimiter = '\t')
##    reader2 = csv.reader(ifile2, delimiter = '\t')

ifiletf = open('/Users/xclair/Documents/google-python-exercises/hg19/TF_information.txt','rt')

readertf = csv.reader(ifiletf, delimiter = '\t')

## generate a list of tuples (motif_id,tf)
mylist1 = list([(row[3],row[6]) for row in readertf]) 


##    ## generate a list of tuples (nm#,gene_name)
##    mylist2 = []
##    mylist2 = list([(row[1],row[-4]) for row in reader2]) 
##    ##print mylist2
##    ##print len(mylist2)
##
##store motif_id, tf in mydict1

##mydict1 = {}
##for k,v in mylist1:
##    if not k in mydict1:
##            mydict1[k] = v
##    else:
##            mydict1[k] += ', '
##            mydict1[k] += v

tfdict = {}
for m,n in mylist1:
    if not n in tfdict:
        tfdict[n]= m
    else:
        tfdict[n] += ','
        tfdict[n] += m
    

##keylist = mydict1.keys()
##keylist.sort()
##d1 = mydict1.items()    

tflist = tfdict.keys()

##
##if __name__=='__main__':
##    sys.exit(main(sys.argv))

ifile = open('/Users/xclair/Documents/google-python-exercises/fimo052912.txt','rt')
reader = csv.reader(ifile, delimiter = '\t')

##motiflist=[]
##for row in reader:
##    if not row[0] in motiflist:
##        motiflist.append(row[0])
##print len(motiflist)

tfname= tfdict.items()
##print tflist

for n in tflist:
    if not tfdict[n] == '.': 
        ofile = open(n + '.gtf','w')
        writer = csv.writer(ofile, delimiter = '\t')
        if ',' in tfdict[n]:
            motifs = tfdict[n].split(',')
            #print motifs
            for m in motifs:    
                ifile.seek(0)
                for row in reader:
                    if row[0] == m:
                        #print m
                        row.append(row[-1])
                        nline = row[1].split('_')
                        if len(nline) == 7:
                            row[-2] = row[3]
                            row[3] = str(int(nline[-2])+int(row[2])-1)#start
                            pval = row[-3]
                            row[-3] = row[4]#strand
                            row[5] = pval#pval
                            row[4] = str(int(nline[-2])+int(row[-2])-1)#end        
                            row[2] = 'NM_' + nline[1]#NM#
                            row[1] = row[0]#motifID
                            row[0] = nline[4]#chr#
                            row[7] = '.'#frame 
                            writer.writerows([row])
                ofile.close()
                
        else:
            ifile.seek(0)
            for row in reader:
                if row[0] == tfdict[n]:
                    #print tfdict[n]
                    row.append(row[-1])
                    nline = row[1].split('_')
                    if len(nline) == 7:
                        row[-2] = row[3]
                        row[3] = str(int(nline[-2])+int(row[2])-1)#start
                        pval = row[-3]
                        row[-3] = row[4]#strand
                        row[5] = pval#pval
                        row[4] = str(int(nline[-2])+int(row[-2])-1)#end        
                        row[2] = 'NM_' + nline[1]#NM#
                        row[1] = row[0]#motifID
                        row[0] = nline[4]#chr#
                        row[7] = '.'#frame 
                        writer.writerows([row])
            ofile.close()



