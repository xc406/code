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

    ##input motif files from the Hughes database
    for infile in glob.glob(os.path.join(path, "M*_0.62.txt")):
        (PATH,FILENAME) = os.path.split(infile)
        (ShortName, Extension)=os.path.splitext(FILENAME)
        print 'Current file being processed is: '+ ShortName

        f = open(infile,'rt')

    ##read file to array
        str1 = f.read()
        list1 = str1.split('\n')
        list0 = []
        ll= len(list1)-1
        if ll > 2:
            for i in range(ll):
                list0.append(list1[i].split('\t'))
            a0 = np.array(list0)
            for i in range(5):
                a0[0,][i]+= ':'
            a1 = np.transpose(a0)
            a2 = np.delete(a1, 0, 0)

    ##convert array to list
            list2 = a2.tolist()

    ##formating list
            list3 = [j for i in list2 for j in i]
    ##  import itertools
    ##  list(itertools.chain(d))
            s1 = '\t'.join(list3[:ll])
            s2 = '\t'.join(list3[ll:(2*ll)])
            s3 = '\t'.join(list3[(2*ll):(3*ll)])
            s4 = '\t'.join(list3[(-1*ll):])
            s = [s1,s2,s3,s4]
            n = '\n'.join(s)

    ##output uniprobe format motif files
            outfile = open(ShortName + '_' +'output' + Extension, 'wt')
            outfile.write(ShortName +'\n'+ n)
            outfile.close()

if __name__=='__main__':
    sys.exit(main(sys.argv))

