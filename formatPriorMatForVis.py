import sys
import os
import csv
import numpy

def main(argv):
    if len(argv) < 2:
        sys.stderr.write("Usage: %s score_matrix \n" % argv[0])
        return 1
    if not os.path.isfile(argv[1]):
        sys.stderr.write('Error: score_matrix %r was not found!\n' % argv[1])
        return 1
    #if not os.path.isfile(argv[3]):
        #sys.stderr.write('Error: nm#_conversion_file %r was not found!\n' % argv[3])
        #return 1

    ##write meme/fimo output to fimoutall by column duplicating 1st and 2nd columns 
    infile = sys.argv[1]

    (path,fname) = os.path.split(infile)

    ifile = open(infile,'rt')
    ofile = open(os.path.join(path, "priorMatVis"),'w')
    reader = csv.reader(ifile, delimiter = '\t')
    writer = csv.writer(ofile, delimiter = '\t')

    #ofile.write('node\tedge\tscore\n')
    ifile.seek(0)
    for column in zip(*reader):
	tfs = column[1:]
	break
    ifile.seek(0)
    for row in reader:
        targets = row
        break
    #print len(targets)#,len(tfs)
    i = 0
    for row in reader:
	for j in xrange(len(targets)):
	    line = [tfs[i],targets[j],row[j+1]]
	    print line
    	    writer.writerows([line])
	i += 1
	    #ofile.write('\n')

    ofile.close()

if __name__=='__main__':
    sys.exit(main(sys.argv))
   
