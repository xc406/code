import sys
import os
import csv
from multiprocessing import Pool, Process

def chunks(tflist,n):
    for i in xrange(0, len(tflist), n):
        yield tflist[i:i+n]##break tflist into chunks

def gffmod(ifile, writer, motif, locidict):
    
    reader = csv.reader(ifile, delimiter = '\t')
    ifile.seek(0)
    for row in reader:
        if row[0] == motif:
            #print row[4]
            #row.append(row[-1])
            #nline = row[1].split('_')
            pstop = row[3]#store stop pos
            row[3] = str(int(locidict[row[1]][1])+int(row[2])-1)#start
            pval = row[-3]#store pval
            row[6] = row[4]#strand
            row[5] = pval#pval
            row[4] = str(int(locidict[row[1]][1])+int(pstop)-1)#stop        
            #nm = nline[0] + '_' + nline[1]#NM#
            motifID = row[0]#motifID
            row[0] = locidict[row[1]][0]#chr
            row[7] = 'gene_id ' + row[1] + '; sequence ' + row[-1]#store group
            row[8] = row[7]#group
            row[1] = motifID
            row[7] = '.'#frame
            row[2] = 'motif'#feature
            writer.writerows([row])
    

def gffgen(tfchunk, tfdict, locidict, ifile):
    
    for n in tfchunk:
        if not tfdict[n] == '.':
        #if n == 'Stat3': #or n == 'IRF4' or n == 'RORC': 
            ofile = open('/home/xc406/data/gff1e3ud10/' + n + '.gff','w')
            writer = csv.writer(ofile, delimiter = '\t')
            if ',' in tfdict[n]:
                motifs = tfdict[n].split(',')
                #print motifs
                for m in motifs:    
                    gffmod(ifile,writer,m,locidict)
                    
            else:
                gffmod(ifile,writer,tfdict[n],locidict)

            ofile.close()
    
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
    
    ##writer = csv.writer(ofile, delimiter = '\t')

    ##(path1,fname1) = os.path.split(infile1)
    ##(path2,fname2) = os.path.split(infile2)

    ifiletf = open('/home/xc406/data/mm9motifs/mm9_motif_output/TF_Information_all_motifs.txt','rt')
    readertf = csv.reader(ifiletf, delimiter = '\t')

    ifileloci = open('/home/xc406/data/mm9genesplus10kb','r')
    readerloci = csv.reader(ifileloci, delimiter = '\t')
    
    mylist1 = list([(row[3],row[6]) for row in readertf])##list of (tf,motifID)
    mylist2 = list([(row[0],int(row[1]),row[-1]) for row in readerloci])##list of (chr,start,gname)
   
    ##generate a dictionary of [tf,motifID]
    tfdict = {}
    for m,n in mylist1:
        if not n in tfdict:
            tfdict[n]= m
        else:
            tfdict[n] += ','
            tfdict[n] += m

    ##generate a dictionary of [gname,(chr,start)]
    locidict = {}
    for c,l,g in mylist2:
	if not g in locidict:
	    locidict[g]= (c,l)
    
    #print len(locidict)    
    tflist = tfdict.keys()##tf list
    gnlist = locidict.keys()##gname list


    tfname= tfdict.items()
    ##print len(tflist) 
    ##print tfdict

    tfchunks = list(chunks(tflist, 500))

    ##print len(tfchunks)
    
    for i in range(len(tfchunks)):
        p = Process(target=gffgen, args=(tfchunks[i],tfdict,locidict,ifile))
        p.start()
        p.join()

if __name__=='__main__':
    #pool = Pool(processor = 12)
    #result = pool.apply_async(f, [])
    sys.exit(main(sys.argv))




