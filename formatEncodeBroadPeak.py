###this script format Encode Chip-Seq output
import sys
import os
import csv
#import operator

def main(argv):
    if len(argv) < 2:
        sys.stderr.write("Usage: %s broadpeak_file\n" % argv[0])
        return 1
    if not os.path.isfile(argv[1]):
        sys.stderr.write('Error: txt_file %r was not found!\n' % argv[1])
        return 1

    #ofile = open('/home/xc406/data/bed_test/' +  shortname + '_o', 'w')
    #writer = csv.writer(ofile, delimiter = '\t')    
    infile = sys.argv[1]

    (path,fname) = os.path.split(infile)
    (shortname, extension) = os.path.splitext(fname)
    ifile = open(infile,'rU')
    reader = csv.reader(ifile, delimiter = '\t')
    
    ofile = open('/home/xc406/data/encode/chipseq/' +  shortname + '_new', 'wt')
    writer = csv.writer(ofile, delimiter = '\t')

    ## macs format [chr,start,end,length,summit,tags,-10*log10(pvalue),fold_enrichment,FDR(%)]
    writer.writerow(['chr','start','end','length','summit','tags','-10*log10(pvalue)','fold_enrichment','FDR(%)'])
    for row in reader:
	row[3],row[4],row[5] = int(row[2])-int(row[1]),(int(row[2])-int(row[1]))/2,row[4]
	writer.writerows([row])    

    ofile.close()

if __name__=='__main__':
    sys.exit(main(sys.argv))



