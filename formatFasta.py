import sys, os

def main(argv):
    if len(argv) < 2:
        sys.stderr.write("Usage: %s <fasta.fa>\nRemove repetitive sequences from fasta files\n" % argv[0])
        return 1
    if not os.path.isfile(argv[1]):
        sys.stderr.write('Error: fasta file %r was not found!\n' % argv[1])
        return 1

    infile = sys.argv[1]
    (path,fname) = os.path.split(infile)
    (shortname,ext) = os.path.splitext(fname) 

    ifile = open (infile,'r')

    f = ifile.readlines()
    l = len(f)
    ##num = l/21

    ##print num
    
    gname = []
    fnew = []
    gnum = 0
    for n in range(l):
	#print f[n]
        if '>' in f[n]:
	    if not 'NNN' in f[n+1]:
                gname.append(f[n])
                gnum+=1
    print gnum

    for i in range(l):
	if f[i] in gname:
	    fnew.append(f[i])
	    for j in range(300):
	        if i+j+1 >= l:
		    break
		elif '>' in f[i+j+1]:
		    break
		else:
		    fnew.append(f[i+j+1])
		    #print f[i+j+1]
	   
    #print fnew

    ##print f

    f2 = ''.join(fnew)

    ofile = open(os.path.join(path,shortname+'norepeat.fa'),'w')
    ofile.write(f2)
    ofile.close


if __name__=='__main__':
    sys.exit(main(sys.argv))


