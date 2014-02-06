### this script parse sam dgf alignment files to leave only 5' cut site
import sys
import os
import csv
#from ftplib import FTP
#import StringIO
import pysam
from os import getenv

#ftp = FTP('hgdownload.cse.ucsc.edu')
#ftp.login()

#sio = StringIO.StringIO()
#def handle_binary(more_data):
#    sio.write(more_data)

#resp = ftp.retrbinary("RETR ", callback = handle_binary)
#sio.seek(0)
#samfile = pysam.Samfile(".bam","rb")

def cutSam(reader,ofile):
    writer = csv.writer(ofile, delimiter = '\t')
    for row in reader:
        if not '@' in row[0]:
            if 'M' in row[5]:
                l = int(row[5][:-1])
                ll = len(row[9])
                if not l == ll:
                    l = ll
                    print 'quality error', row
                lq = len(row[10])
                if not ll == lq:
                    print 'sequence quality inconsistent', row
                row[5] = '1M'
                if not row[10][0] == '*':
                    if len(row[10][0]) == 1:
                        row[10] = row[10][0]
                    else:
                        row[10] = '0'
                else:
                    #print 'quality score missing', row[0]
                    row[10] = '0'
            ##reads mapped to reverse strand
                if row[1] == '16':
                    row[9] = row[9][l-1]
                    row[3] = str(int(row[3]) + l - 1)
            ##reads mapped to forward strand
                if row[1] == '0':
                    row[9] = row[9][0]
            #else:
                #print 'read length == ', row[5]
            #print row
        #else:
            #print row[0]
        writer.writerows([row])
    ofile.close()

def cutBam(bamfile,ofile):
    for s in bamfile:
        a = pysam.AlignedRead()
        a.qname = s.qname
        a.flag = s.flag
        a.rname = s.rname
        a.mapq = s.mapq
        a.tags = s.tags
        a.mrnm = s.mrnm
        a.mpos = s.mpos
        a.is_paired = False
        a.cigar = [(0,1)]
        if s.flag == 16:
            a.pos = int(s.pos)+int(s.rlen)-1
            a.seq = s.seq[s.rlen-1]
	    if not s.qual[s.rlen-1] == '*':
                a.qual = s.qual[s.rlen-1]
		#if duplicateCounter[(a.rname,a.pos)] < cut:
		ofile.write(a)
	    #else:
		#a.qual = '0'
        elif s.flag == 0:
            a.pos = s.pos
            a.seq = s.seq[0]
	    if not s.qual[0] == '*':
                a.qual = s.qual[0]
	    #else:
		#a.qual = '0'
		#if duplicateCounter[(a.rname,a.pos)] < cut:
                ofile.write(a)
    ofile.close()

def cleanBam(bamfile,ofile,duplicateCounter,cut):
    for s in bamfile:
	if duplicateCounter[(s.rname,s.pos)] < cut:
	    ofile.write(a)
    ofile.close()

def countDuplicate(bamfile):
    duplicateCounter = defaultdict(int)
    for s in bamfile:
        #print s.rname
        #print s.pos
        #if not (s.rname,s.pos) in duplicateCounter:
            #print (s.rname,s.pos)
        duplicateCounter[(s.rname,s.pos)] += 1
    return duplicateCounter

def main(argv):
    if len(argv) < 2:
        sys.stderr.write("Usage: %s bam_file\n" % argv[0])
        return 1
    if not os.path.isfile(argv[1]):
        sys.stderr.write('Error: bam_file %r was not found!\n' % argv[1])
        return 1

    if getenv('MYTMP'):
        opath = getenv('MYTMP')
        #print opath, 'output path'

    #ofile = open('/home/xc406/data/bed_test/' +  shortname + '_o', 'w')
    #writer = csv.writer(ofile, delimiter = '\t')    
    infile = sys.argv[1]

    (path,fname) = os.path.split(infile)
    (shortname, extension) = os.path.splitext(fname)

    #ifile = open(infile,'rt')
    #reader = csv.reader(ifile, delimiter = '\t')
    
    #ofile = open(os.path.join('/data/cgsb/bonneau/xchen/dhs/', shortname + 'Cut5.sam'), 'wt')
    #writer = csv.writer(ofile, delimiter = '\t')

    bamfile = pysam.Samfile(infile,'rb')
    ofile = pysam.Samfile(os.path.join(opath,shortname+'Cap.sam'),'wh', template = bamfile)
    duplicateCounter = countDuplicate(bamfile)
    cleanBam(bamfile,ofile,duplicateCounter,150)

if __name__=='__main__':
    sys.exit(main(sys.argv))



