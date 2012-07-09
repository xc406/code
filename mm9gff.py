import sys
import os
import csv

def main(argv):
    if len(argv) < 2:
        sys.stderr.write("Usage: %s fimo_output_file TF_Info_file nm#_conversion_file\n" % argv[0])
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

    ifiletf = open('/home/xc406/data/motif_output/TF_Information_all_motifs.txt','rt')
    readertf = csv.reader(ifiletf, delimiter = '\t')

    mylist1 = list([(row[3],row[6]) for row in readertf])

    tfdict = {}
    for m,n in mylist1:
        if not n in tfdict:
            tfdict[n]= m
        else:
            tfdict[n] += ', '
            tfdict[n] += m

    tflist = tfdict.keys()

##ifile = open('/home/xc406/data/fimo052912.txt','rt')
##reader = csv.reader(ifile, delimiter = '\t')

    tfname= tfdict.items()
    ##print tflist

    for n in tflist:
        if not tfdict[n] == '.': 
            ofile = open('/home/xc406/data/mm9gff/' + n + '.gff','w')
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
                                row[7] = 'gene_id ' + row[2] + '; sequence ' + row[-1]#group
                                row[-1] = row[7]#group
                                row[7] = '.'#frame
                                row[2] = 'motif'#feature
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
                            row[7] = 'gene_id ' + row[2] + '; sequence ' + row[-1]#group
                            row[-1] = row[7]#group
                            row[7] = '.'#frame
                            row[2] = 'motif'#feature 
                            writer.writerows([row])
                ofile.close()


if __name__=='__main__':
    sys.exit(main(sys.argv))




