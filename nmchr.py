import sys, os

def main(argv):
    if len(argv) < 2:
        sys.stderr.write("Usage: %s <upstream1000.fa>\n" % argv[0])
        return 1
    if not os.path.isfile(argv[1]):
        sys.stderr.write('Error: fasta file %r was not found!\n' % argv[1])
        return 1

    infile = sys.argv[1]
    (path,fname) = os.path.split(infile)

    ifile = open (infile,'r')

    f = ifile.readlines()
    l = len(f)
    ##num = l/21

    ##print num

    gname = []
    gnum = 0
    for n in range(l):
        if n % 21 == 0:
            gname.append(f[n])
            gnum+=1
    print gname

    for i in range(gnum):
        sgname = gname[i].split(' ')
        gname[i] = sgname[0]
        ##nm = gname[i].split('_up_1000_')
        ##gname[i] = nm[0]
        for n in range((i*21),(i*21+1)):
            if n == i*21:
                f[n] = gname[i] + '\n'
    ##print f
    f2 = ''.join(f)

    ofile = open(os.path.join(path,'nmchrupstream1000.fa'),'w')
    ofile.write(f2)
    ofile.close


if __name__=='__main__':
    sys.exit(main(sys.argv))



