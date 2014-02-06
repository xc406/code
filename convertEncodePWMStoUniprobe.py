import csv
import os
import sys
import re
import glob
import numpy as np

def main(argv):
    if len(argv) < 2:
        sys.stderr.write("Usage: %s Path-to-motif-files\n" % argv[0])
        return 1
    if not os.path.exists(argv[1]):
        sys.stderr.write('Error: Path %r was not found!\n' % argv[1])
	return 1
    
    path = sys.argv[1]

    ##input motif files from ENCODE
    for infile in glob.glob(os.path.join(path, "*_*")):
        (PATH,FILENAME) = os.path.split(infile)
        #(ShortName, Extension)=os.path.splitext(FILENAME)
        print 'Current file being processed is: '+ FILENAME

        f = open(infile,'rt')

    ##read file to array
        str1 = f.read()
        list1 = str1.split('\n')
	list2 = []
	for i in list1:
	    i = i.rstrip()
	    list2.append(i.split(' ')[1:])
	list0 = list2[1:len(list2)-1]
        a0 = np.array(list0)
	a1 = np.transpose(a0)
	a2 = np.flipud(a1)
	a3 = np.fliplr(a2)##flip the matrix to represent the binding site on the opposite strand
	#print a0,a1
	#print a2,a3
	list3 = a3.tolist()
	ntlist = ['A:','C:','G:','T:']
	list4 = []
	f.seek(0)
	info = f.readline().rstrip().rstrip()
	if '\t' in info:
	    tfname = info.split('\t')[1].split('_')[0].upper()
	else:
	    tfname = info.split(' ')[1].split('_')[0].upper()
	print info
	print tfname
	if '/' in tfname:
	    tfname = tfname.replace('/','')
	list4.append(info)
	for i in xrange(len(list3)):
	    row = ntlist[i] + '\t' + '\t'.join(list3[i])
	    list4.append(row)
	s = '\n'.join(list4)
	outfile = open(os.path.join(path + 'test', FILENAME + '_' + tfname + '_' +'uniprobe'), 'wt')
	outfile.write(s)
	outfile.close()
    ##flatten a list of a list
        #list4 = [j for i in list3 for j in i]

    ##  list(itertools.chain(d))
            #s1 = '\t'.join(list3[:ll])
            #s2 = '\t'.join(list3[ll:(2*ll)])
            #s3 = '\t'.join(list3[(2*ll):(3*ll)])
            #s4 = '\t'.join(list3[(-1*ll):])
            #s = [s1,s2,s3,s4]
            #n = '\n'.join(s)

    ##output uniprobe format motif files
            #outfile = open(os.path.join(PATH + '/motifoutput', ShortName + '_' +'uniprobe'), 'wt')
            #outfile.write(ShortName +'\n'+ n)
            #outfile.close()

if __name__=='__main__':
    sys.exit(main(sys.argv))

