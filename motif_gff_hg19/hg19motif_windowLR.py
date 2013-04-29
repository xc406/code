import sys
import os
import csv
#import operator

def main(argv):
    if len(argv) < 2:
        sys.stderr.write("Usage: %s gff_file\n" % argv[0])
        return 1
    if not os.path.isfile(argv[1]):
        sys.stderr.write('Error: gff_file %r was not found!\n' % argv[1])
        return 1
##    if not os.path.isfile(argv[2]):
##        sys.stderr.write('Error: TF_Info_file %r was not found!\n' % argv[2])
##        return 1
##    if not os.path.isfile(argv[3]):
##        sys.stderr.write('Error: nm#_conversion_file %r was not found!\n' % argv[3])
##        return 1
    
    infile = sys.argv[1]

    (path,fname) = os.path.split(infile)
    (tfname,ext) = os.path.splitext(fname)

    ifile = open(infile,'rt')
    reader = csv.reader(ifile, delimiter = '\t')

##ifile = open('/home/xc406/data/fimo052912.txt','rt')
##reader = csv.reader(ifile, delimiter = '\t')

    ofileL = open('/home/xc406/data/hg19gff1e3ud5_final/' + tfname+'_winL10.gff', 'w')
    writerL = csv.writer(ofileL, delimiter = '\t')
    ofileR = open('/home/xc406/data/hg19gff1e3ud5_final/' + tfname+'_winR10.gff', 'w')
    writerR = csv.writer(ofileR, delimiter = '\t')
#		sortedlist = sorted(reader, key=operator.itemgetter(3), reverse=False)	    

    ##expand the motif window to upstream and downstream 100bp
    ##L count
    ifile.seek(0)
    #pline = ''
    for row in reader:
        row[3],row[4] = str(int(row[3])-10),row[3]
	#row[4] = row[3]#str(int(row[4])+100)
	#group = row[-1].split(';')
	#gid = group[0].split(' ')
	#for k in nmlist:
	    #if gid[1] == k:
		#gid[1] = gndict[k]
		#group[0] = gid[0] + ' ' + gid[1]
		#group = group[0] + ';' + group[1]
  		#row[-1] = group
		#if not row == pline: 
	writerL.writerows([row])
    ofileL.close()

    ##R count
    ifile.seek(0)
    for row in reader:
	#row[3] = row[4]
	row[3],row[4] = row[4],str(int(row[4])+10)
	writerR.writerows([row])
                    #pline = row
    ofileR.close()


if __name__=='__main__':
    sys.exit(main(sys.argv))



