import sys, os, csv
import numpy as np
import scipy
from collections import defaultdict
import operator
import pickle

def matchem(infile, emdict):
    ifile = open(infile, "rt")
    reader = csv.reader(ifile, delimiter = '\t')
    ifile.seek(0)
    for row in reader:
	if '_' in row[1]: 
	    #print row
	    seq_len = int(row[1].split('_')[2]) - int(row[1].split('_')[1])
	    motif_len = int(row[3])-int(row[2])+1
	    emdict[(row[1],row[0].split('_')[0])] += 1.0/seq_len
    return emdict  
    
def JaccardIndex(infile, mdict):
    ifile = open(infile, "rt")
    reader = csv.reader(ifile, delimiter = '\t')
    ifile.seek(0)
    for row in reader:
	mo = row[0].split('_')[0]
	e = row[1]
	if 'M' in mo:
	    if not e in mdict[mo]:
	        mdict[mo].append(e)
    ji = defaultdict(float)
    for m1 in mdict:
	for m2 in mdict:
	    ji[(m1,m2)] = len(intersect(m1,m2))/float(len(union(m1,m2)))
    return ji	    

def intersect(a, b):
     """ return the intersection of two lists """
     return list(set(a) & set(b))

def union(a, b):
    """ return the union of two lists """
    return list(set(a) | set(b))

def flatten(lst):
    for elem in lst:
        if type(elem) == list:
            for item in elem:
                yield str(item)
        else:
            yield str(elem)

def main(argv):
    if len(argv) < 2:
        sys.stderr.write("Usage: %s fimo_output_file \n" % argv[0])
        return 1
    if not os.path.isfile(argv[1]):
        sys.stderr.write('Error: fimo_output_file %r was not found!\n' % argv[1])
        return 1

    infile = sys.argv[1]

    (path,fname) = os.path.split(infile)

    ofile = open(os.path.join(path,'JIMatrix'),'wt')
    writer = csv.writer(ofile, delimiter = '\t')

    ##(path1,fname1) = os.path.split(infile1)
    ##(path2,fname2) = os.path.split(infile2)

    ifiletf = open('/home/xc406/data/mm9motifs80/TF_Information80mm9.txt','rt')
    readertf = csv.reader(ifiletf, delimiter = '\t')
    
    ifilemotif = open('/home/xc406/data/mm9motifs80/mappable/motifoutput/mm9motifs80.txt', 'rt')
    readermotif = csv.reader(ifilemotif,delimiter = '\t')
    	
    mylist1 = list([(row[3],row[6]) for row in readertf])##list of (motifID,tf)

    ##generate a dictionary of [tf,motifID]
    emdict = defaultdict(float)
    mdict = defaultdict(list)
    datadict = defaultdict(list)
    tfdict = {}
    for m,n in mylist1:
        if not n in tfdict:
            tfdict[n]= m
        else:
            tfdict[n] += ','
            tfdict[n] += m
    
    tflist = tfdict.keys()##tf list

    tfname= tfdict.items()
 
    jidict = JaccardIndex(infile,mdict)

    #print jidict

    with open('motif2tfname.txt','r') as f:
        mydict1 = pickle.load(f)
    with open('emodict.txt','r') as f:
        emodict = pickle.load(f)
    print emodict
    ##create motif/element list
    mlist = []
    for row in readermotif:
	motifid = row[0].split(' ')[-1].split('_')[0]
	if not motifid in mlist:
	    temp = [motifid,mydict1[motifid]]
	    mlist.append(('_').join(temp))
    
    mmlist = []
    #print len(emdict_final)
    for k,v in jidict:
	if not k in mmlist:
	    mmlist.append(k)
    #print len(elist)
	
    nmlist = []
    for mn1 in mlist:
	m1 = mn1.split('_')[0]
	n1 = mn1.split('_')[1]
	for mn2 in mlist:
	    m2 = mn2.split('_')[0]
	    n2 = mn2.split('_')[1]
	    if n1 in emodict:
		if n2 in emodict:
		    if (m1 + '_0.80' in emodict[n1]) and (m2 + '_0.80' in emodict[n2]):
		        if not [m1,n1] in nmlist:
		            nmlist.append([m1,n1])
	                if (m1,m2) in jidict:
		            datadict[(m1,n1)].append((m2,jidict[(m1,m2)]))
	                else:
		            datadict[(m1,n1)].append((m2,0.0))

    nmlist.sort(key=operator.itemgetter(0))
    nlist = []
    for nm in nmlist:
	nlist.append(nm[0] + '_' + nm[1])
    #print nlist
    #print list(flatten(['seq',mlist]))
    writer.writerows([list(flatten(['seq', nlist]))])
	
    rownamelist = []
    for (m1,n1) in datadict:
	rownamelist.append((m1,n1))	        
        datadict[(m1,n1)].sort(key=operator.itemgetter(0))
    #print datadict	
    rownamelist.sort(key=operator.itemgetter(0))
    print rownamelist

    for (m1,n1) in rownamelist:
	jilistsort = []
	jilistsort.append(m1 + '_' + n1)
	for m2,ji in datadict[(m1,n1)]:
	    jilistsort.append(ji)
	#print len(nlistsort)
        writer.writerows([jilistsort])

    #print datadict
    #emmatrix = emdict_final#np.ndarray(shape=(3789,651), dtype=float, order='C')
    #emmatrix.fill(0.0)
    #print emmatrix
    ofile.close()

if __name__=='__main__':
    sys.exit(main(sys.argv))




