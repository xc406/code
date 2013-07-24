import sys
import os
import csv

def main(argv):
    if len(argv) < 2:
        sys.stderr.write("Usage: %s TF_Info_file\n" % argv[0])
        return 1
    if not os.path.isfile(argv[1]):
        sys.stderr.write('Error: TF_Info_file %r was not found!\n' % argv[1])
        return 1
##    if not os.path.isfile(argv[2]):
##        sys.stderr.write('Error: TF_Info_file %r was not found!\n' % argv[2])
##        return 1
##    if not os.path.isfile(argv[3]):
##        sys.stderr.write('Error: nm#_conversion_file %r was not found!\n' % argv[3])
##        return 1
    
    infile = sys.argv[1]

    (path,fname) = os.path.split(infile)

    ifile = open(infile,'rt')
    reader = csv.reader(ifile, delimiter = '\t')
    ##writer = csv.writer(ofile, delimiter = '\t')

    ##(path1,fname1) = os.path.split(infile1)

    #ifiletf = open('/home/xc406/data/Homo_sapiens_2013_02_20_3-54_pm/TF_Information_all_motifs71.txt','rt')
    #readertf = csv.reader(ifiletf, delimiter = '\t')

    #(path,fname) = os.path.split(ifiletf)
    
    mtlist = list([(row[3],row[6],row[15]) for row in reader])

    tfdict = {}
    for m,t,mt in mtlist:
	#if (mt == 'PBM') or (mt == 'SELEX'):
	if mt == 'Jolma':
            if not t in tfdict:
                tfdict[t]= m
            else:
                tfdict[t] += ','
                tfdict[t] += m

    tflist = tfdict.keys()
    #tflist = ["GATA3","RORC","TBX21","FOXP3","BATF","IRF4","IRF8"]
    tfname= tfdict.items()
    ##print tflist 
    ##print tfdict
  
    for t in tflist:
        motiflist = []
	if not tfdict[t] == '.':
	#if n == 'RORC': #or n == 'IRF4' or n == 'RORC': 
            if ',' in tfdict[t]:
                motifs = tfdict[t].split(',')
                for m in motifs:
		    motiflist.append(m)
	    else:
		motiflist.append(tfdict[t])
            #print motiflist
    	    outfile = open(os.path.join(path,'HTSELEX/'+t+'.txt'), 'wt')
            for m in motiflist:
	        outfile.write(m +'\n')
            outfile.close()

if __name__=='__main__':                    
    sys.exit(main(sys.argv))




