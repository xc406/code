import sys
import csv
import os

def main(argv):
    if len(argv) < 2:
        sys.stderr.write("Usage: %s gff_file \n" % argv[0])
        return 1
    if not os.path.isfile(argv[1]):
        sys.stderr.write('Error: gff_file %r was not found! \n' % argv[1])
        return 1

    infile = sys.argv[1]
    
    (path,fname) = os.path.split(infile)

    ifile = open(infile,'rt')
    ofile = open(os.path.join(path, '/mm9gff_update/' + fname),'w')
    reader = csv.reader(ifile, delimiter = '\t')
    writer = csv.writer(ofile, delimiter = '\t')

#inputfile = open('/Users/xichen/Downloads/htseq_output1', 'r')
#readerfile = csv.reader(inputfile, delimiter = '\t')

#inputfile1 = open('/home/xc406/Downloads/Arx.gff', 'rt')
#readerfile1 = csv.reader(inputfile1, delimiter = '\t')
#ofile = open( 'Arx_update.gff' , 'w' )
#writerfile1 = csv.writer(ofile, delimiter = '\t')
    previous = ''
    previous_start = 0
    previous_stop = 0
    previous_chr = ''
    previous_strand = ''
    counter = 0
    geneidlist=[]
    for row in ifile:
	if row[0] == previous_chr & row[6] == previous_strand:
	    if int(row[3]) in [previous_start, previous_stop]:
	    		


	 
	line = row.split('fimo\t')  
        list1 = line[-1].split(' ')
        geneid = list1[1][0:-1]
        if geneid == previous:
            counter +=1
            geneid = geneid + '_' + str(counter)
        else:
            previous = geneid
            counter = 0
    #    print geneid    
        list1[1] = geneid + ';'
        list1[3] = list1[3][0:-2]
        line = line[0:-1]
        line.append(list1[0])
        line.append(list1[1])
        line.append(list1[2])
    #    print list1[3]
    #    print line
        line.append(list1[3])
        #print line
    #    line.append('\n')
        writer.writerows([line])

    ofile.close()

if __name__=='__main__':
    sys.exit(main(sys.argv))

#nmlist = []
#nmlistnormalized = []
#for row in readerfile:
#    row[1] = int(row[1])
##    if 'NM_' in row[0]:
##        if not row[1] == 0:
##            nmlist.append(row[1])
#    tot_count += row[1]

#print tot_count
#inputfile.seek(0)
#for row in readerfile:
#    row[1] = int(row[1])  
#    if 'NM_' in row[0]:
#        if not row[1] == 0:
#            nmlistnormalized.append(float(row[1])/tot_count)

##print nmlistnormalized
