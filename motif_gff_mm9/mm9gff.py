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

    ifile = open(infile,'rt')
    reader = csv.reader(ifile, delimiter = '\t')
    ##writer = csv.writer(ofile, delimiter = '\t')

    ##(path1,fname1) = os.path.split(infile1)
    ##(path2,fname2) = os.path.split(infile2)

    ifiletf = open('/home/xc406/data/mm9motifs/mm9_motif_output/TF_Information_all_motifs.txt','rt')
    readertf = csv.reader(ifiletf, delimiter = '\t')

    mylist1 = list([(row[3],row[6]) for row in readertf])

    tfdict = {}
    for m,n in mylist1:
        if not n in tfdict:
            tfdict[n]= m
        else:
            tfdict[n] += ','
            tfdict[n] += m

    tflist = tfdict.keys()

##ifile = open('/home/xc406/data/fimo052912.txt','rt')
##reader = csv.reader(ifile, delimiter = '\t')

    tfname= tfdict.items()
    ##print tflist 
    ##print tfdict

    for n in tflist:
        if not tfdict[n] == '.':
	#if n == 'Stat3': #or n == 'IRF4' or n == 'RORC': 
            ofile = open('/home/xc406/data/mm9gff1e3/' + n + '.gff','w')
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
                                pstop = row[3]#store stop pos
                                row[3] = str(int(nline[-2])+int(row[2])-1)#start
                                pval = row[-3]#store pval
                                row[-3] = row[4]#strand
                                row[5] = pval#pval
                                row[4] = str(int(nline[-2])+int(pstop)-1)#stop        
                                nm = nline[0] + '_' + nline[1]#NM#
                                row[1] = row[0]#motifID
                                row[0] = nline[4]#chr# 
                                row[7] = 'gene_id ' + nm + '; sequence ' + row[-1]#store group
                                row[-1] = row[7]#group
                                row[7] = '.'#frame
                                row[2] = 'motif'#feature
                                writer.writerows([row])
                    
            else:
                ifile.seek(0)
                for row in reader:
                    if row[0] == tfdict[n]:
                        #print tfdict[n]
                        row.append(row[-1])
                        nline = row[1].split('_')
                        if len(nline) == 7:
                            pstop = row[3]
                            row[3] = str(int(nline[-2])+int(row[2])-1)#start
                            pval = row[-3]
                            row[-3] = row[4]#strand
                            row[5] = pval#pval
                            row[4] = str(int(nline[-2])+int(pstop)-1)#end        
                            nm = nline[0] + '_' + nline[1]#NM#
                            row[1] = row[0]#motifID
                            row[0] = nline[4]#chr#
                            row[7] = 'gene_id ' + nm + '; sequence ' + row[-1]#group
                            row[-1] = row[7]#group
                            row[7] = '.'#frame
                            row[2] = 'motif'#feature 
                            writer.writerows([row])
            ofile.close()


if __name__=='__main__':
    sys.exit(main(sys.argv))




