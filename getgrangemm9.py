import sys
import os
import csv
#import operator

def main(argv):
    if len(argv) < 3:
        sys.stderr.write("Usage: %s Refseq_file Fasta_file\n" % argv[0])
        return 1
    if not os.path.isfile(argv[1]):
        sys.stderr.write('Error: Refseq_file %r was not found!\n' % argv[1])
	return 1
    if not os.path.isfile(argv[2]):
        sys.stderr.write('Error: Fasta_file %r was not found!\n' % argv[2])
        return 1
    
    infile = sys.argv[1]
    infile1 = sys.argv[2]

    (path,fname) = os.path.split(infile)
    (shortname, extension) = os.path.splitext(fname)
    ifile = open(infile,'rU')
    ifile1 = open(infile1,'rU')
 
    reader = csv.reader(ifile, delimiter = '\t')
    reader1 = csv.reader(ifile1, delimiter = '\t')
    ##writer = csv.writer(ofile, delimiter = '\t')

    ##(path1,fname1) = os.path.split(infile1)
    ##(path2,fname2) = os.path.split(infile2)

##ifile = open('/home/xc406/data/fimo052912.txt','rt')
##reader = csv.reader(ifile, delimiter = '\t')

    ofile = open('/home/xc406/data/' + 'mm9ud5.bed', 'wt')
    writer = csv.writer(ofile, delimiter = '\t')
	
#sortedlist = sorted(reader, key=operator.itemgetter(3), reverse=False)	    
    
    chromlist = []
    chromdict = {}
    for row in reader1:
        if '>' in row[0]:
            i = row[0][1:]
            if not i in chromlist:
                chromlist.append(i)
                chromdict[i] = 0
        if not '>' in row[0]:
            chromdict[i] += len(row[0]) ## create a dictionary of chr and size
    
    ##write row_number/target(mtl.id) to a list
    genelist = []

    ifile.seek(0)
 
    for row in reader:
	chrom = row[2]
	#if not row == pline:
	if chrom in ['chrM','chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chrX','chrY']:
	    gname = row[-4]
	    if not gname in genelist:
	        genelist.append(gname)
	        newrow = ['','','','']
	        newrow[0] = chrom##chromosome name
		strand = row[3]
		tss = row[4]
		tes = row[5] 
		#if (int(row[5])+10000<=chromdict[row[2]]) and (int(row[4])-10000>=1):##upstream and downstream 10kb
	        if strand == '+':
		    if int(tss)-5000>=1:
		        newrow[1] = str(int(tss)-5000)#start
		    else:
		        newrow[1] = '1'
         	    if int(tss)+4999<=chromdict[chrom]:##TSS upstream 10kb and downstream 1kb 
	                newrow[2] = str(int(tss)+4999)#end
	            else:
		        newrow[2] = str(chromdict[chrom])##make sure 10kb upstream or downstream don't go out of chromosomal range
		else:
		    if int(tes)-5000>=1:
                        newrow[1] = str(int(tes)-5000)#start
                    else:
                        newrow[1] = '1'
                    if int(tes)+4999<=chromdict[chrom]:##TSS upstream 10kb and downstream 1kb
                        newrow[2] = str(int(tes)+4999)#end
                    else:
                        newrow[2] = str(chromdict[chrom])##make sure 10kb upstream or downstream don't go out of chromosomal range
		newrow[3] = gname##gene name  
                writer.writerows([newrow])

    ofile.close()

if __name__=='__main__':
    sys.exit(main(sys.argv))



