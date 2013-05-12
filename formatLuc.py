###this script format luciferase assay gold standard
import sys
import os
import csv
#import operator

def main(argv):
    if len(argv) < 2:
        sys.stderr.write("Usage: %s txt_file\n" % argv[0])
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
    
    ofile = open('/home/xc406/data/gold_standard/' +  shortname + '_new.bed', 'wt')
    writer = csv.writer(ofile, delimiter = '\t')

    #line = ['','',0,0,'',0.0,0.0]## store gs as [target_gene,chr,start,end,strand,ratio,fdr]
    line = ['',0,0,'',0.0,'']## store gs as [chr,start,end,target_gene,ratio,strand]
    for row in reader:
	line = [row[5].split(':')[0],row[5].split(':')[1].split('-')[0],row[5].split(':')[1].split('-')[1],row[1],row[-2],row[-3]]
	writer.writerows([line])    

    ofile.close()

if __name__=='__main__':
    sys.exit(main(sys.argv))



