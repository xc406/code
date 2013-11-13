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
	    #seq_len = int(row[1].split('_')[2]) - int(row[1].split('_')[1])
	    motif_len = int(row[3])-int(row[2])+1
	    em_id = (row[1],row[0].split('_')[0])##id in the form of seq_motif
	    #if not em_id in emdict:
	    emdict[em_id] += 1.0#float(row[6])#1.0#/seq_len
	    #else:
		#emdict[em_id] *= float(row[6])##taking the product of p-vals
    return emdict  
    
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

    ofile = open(os.path.join(path,'fireMatrixECleanExpr'),'wt')
    writer = csv.writer(ofile, delimiter = '\t')

    ##(path1,fname1) = os.path.split(infile1)
    ##(path2,fname2) = os.path.split(infile2)

    ifiletf = open('/home/xc406/data/mm9motifs80/TF_Information80mm9.txt','rt')
    readertf = csv.reader(ifiletf, delimiter = '\t')
    
    ifilemotif = open('/home/xc406/data/mm9motifs80/mappable/motifoutput/mm9motifs80.txt', 'rt')
    readermotif = csv.reader(ifilemotif,delimiter = '\t')
    	
    mylist1 = list([(row[3],row[6]) for row in readertf])##list of (motifID,tf)

    """list of enriched tfs filtered by expression in esc """
    exprlist = ['Klf4','Sox2','Pou5f1','Sp1','Zfx','Zic3','Rfx2','Zfp281','Tcfap2c','Zic2','Tcf3','Zbtb7b']
    
    ##generate a dictionary of [tf,motifID]
    emdict = defaultdict(float)
    datadict = defaultdict(lambda : defaultdict(float))
    tfdict = {}
    for m,n in mylist1:
        if not n in tfdict:
            tfdict[n]= m
        else:
            tfdict[n] += ','
            tfdict[n] += m
    
    tflist = tfdict.keys()##tf list

    tfname= tfdict.items()
 
    emdict_final=matchem(infile,emdict)
    #print emdict_final

    with open('motif2tfname.txt','r') as f:
        mydict1 = pickle.load(f)
    with open('emodict.txt','r') as f:
        emodictp = pickle.load(f)
    
    emodict = defaultdict(list)

    for k in emodictp:
        if k in exprlist:
           emodict[k] = emodictp[k]

    print emodict
    ##create motif/element list
    mlist = []
    for row in readermotif:
	motifid = row[0].split(' ')[-1].split('_')[0]
	if not motifid in mlist:
	    temp = [motifid,mydict1[motifid]]
	    mlist.append(('_').join(temp))
    
    elist = []
    #print len(emdict_final)
    for k,v in emdict_final:
	if not k in elist:
	    elist.append(k)
    #print len(elist)
    nlist = []
    for e in elist:
	for mn in mlist:
	    m = mn.split('_')[0]
	    n = mn.split('_')[1]
	    if n in emodict:
		if m + '_0.80' in emodict[n]:
		    if not n in nlist:
		        nlist.append(n)
	            if (e,m) in emdict_final:
			if not e in datadict:
		            datadict[e][n] = emdict_final[(e,m)]
			elif not n in datadict[e]:
			    datadict[e][n] = emdict_final[(e,m)]
			else:
			    datadict[e][n] += emdict_final[(e,m)]
	            else:
		        datadict[e][n] = 0.0

    #nmlist.sort(key=operator.itemgetter(1))
    nlist.sort()
    #for nm in nmlist:
	#nlist.append(nm[1])
    print nlist
    #print list(flatten(['seq',mlist]))
    writer.writerows([list(flatten(['seq', nlist]))])
    
    rownamelist = []
    for e in datadict:
	rownamelist.append(e)
        datadict[e] = sorted(datadict[e].items())	        
        #datadict[e].sort(key=operator.itemgetter(0))
    
    rownamelist.sort()
	
    for e in rownamelist:
	nlistsort = []
	nlistsort.append(e)
	for (n,count) in datadict[e]:
	    nlistsort.append(float(count)/len(emodict[n]))
	print nlistsort
	if not sum(nlistsort[1:]) == 0.0:##discard low signal sequences
	    writer.writerows([nlistsort])

    #print datadict
    #emmatrix = emdict_final#np.ndarray(shape=(3789,651), dtype=float, order='C')
    #emmatrix.fill(0.0)
    #print emmatrix
    ofile.close()

if __name__=='__main__':
    sys.exit(main(sys.argv))




