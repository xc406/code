import csv
import os
import sys
import re
import glob
import numpy as np

##!! missing dual orientation comparison
def main(argv):
    if len(argv) < 3:
        sys.stderr.write("Usage: %s Path-to-Database-1 Path-to-Database-2\n" % argv[0])
        return 1
    if not os.path.exists(argv[1]):
        sys.stderr.write('Error: Path-1 %r was not found!\n' % argv[1])
	return 1
    if not os.path.exists(argv[2]):
	sys.stderr.write('Error: Path-2 %r was not found!\n' % argv[2])
    	return 1

    path1 = sys.argv[1]
    path2 = sys.argv[2]
    count = 0
    ##input motif files from ENCODE
    for infile1 in glob.glob(os.path.join(path1, "*_uniprobe")):
        (PATH1,FILENAME1) = os.path.split(infile1)
	for infile2 in glob.glob(os.path.join(path2, "*_uniprobe")):
            (PATH2, FILENAME2) = os.path.split(infile2)
            #print 'Comparing: '+ FILENAME1 + ',' + FILENAME2

            f1 = open(infile1,'rt')
	    str1 = f1.read()
	    list1 = str1.split('\n')
	    f2 = open(infile2,'rt')
            str2 = f2.read()
            list2 = str2.split('\n')
	    #print list1[1:]
	    #print list2[1:]
	    if list1[1:] == list2[1:]:
	    	print 'identical' + FILENAME1 + ',' + FILENAME2
		print list1[0]
		os.remove(infile2)
		count +=1
    print count
	#list4 = []
	#f.seek(0)
	#list4.append(f.readline().rstrip().rstrip())
	#for i in xrange(len(list3)):
	#    row = ntlist[i] + '\t' + '\t'.join(list3[i])
	#    list4.append(row)
	#s = '\n'.join(list4)
	#outfile = open(os.path.join(path + '/motifoutput', FILENAME + '_' +'uniprobe'), 'wt')
	#outfile.write(s)
	#outfile.close()
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

