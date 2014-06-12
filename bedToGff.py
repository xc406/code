import sys
import os
import csv
import random

"""been changed around!"""
def main(argv):
    if len(argv) < 2:
        sys.stderr.write("Usage: %s bed_file\n" % argv[0])
        return 1
    if not os.path.isfile(argv[1]):
        sys.stderr.write('Error: bed_file %r was not found!\n' % argv[1])
        return 1
##    if not os.path.isfile(argv[2]):
##        sys.stderr.write('Error: TF_Info_file %r was not found!\n' % argv[2])
##        return 1
##    if not os.path.isfile(argv[3]):
##        sys.stderr.write('Error: nm#_conversion_file %r was not found!\n' % argv[3])
##        return 1
    
    infile = sys.argv[1]

    #(path,fname) = os.path.split(infile)

    ifile = open(infile,'rt')
    reader = csv.reader(ifile, delimiter = '\t')
    ##writer = csv.writer(ofile, delimiter = '\t')

    (path,fname) = os.path.split(infile)
    (name,ext) = os.path.splitext(fname)

    #ifile = open('/home/xc406/data/mm9upstream10kb.bed','rt')
    #reader = csv.reader(ifile, delimiter = '\t')

    ofile = open(os.path.join(path,name+'.gff'),'w')#'/home/xc406/data/hg19_gname_dist.gff', 'w')
    writer = csv.writer(ofile, delimiter = '\t')
    
    line = ['','.','feature',0,0,'1','+','.','']
    for row in reader:
	line[0] = row[0]
	##change from 0-based bed format to 1-based gff format
	line[3],line[4] = int(row[1])+1,int(row[2])+1
	line[-1] = 'gene_id '+ str(row[0])+'_'+str(int(row[1])+1)+'_'+str(int(row[2])+1)#+'_'+str(row[3])+'_'+str(row[4])
	#row[1],row[2],row[3],row[4],row[5] = row[3].split('_')[3],'fire',row[1],row[2],'1'
	#row[1],row[2],row[3],row[4],row[5],row[6],row[7],row[8] = str(row[0]+'_'+row[1]+'_'+row[2]),'TSS',row[1],row[2],'1','+','.',row[6]	
	writer.writerows([line])
    
    ofile.close()

if __name__=='__main__':
    sys.exit(main(sys.argv))




