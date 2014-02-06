import csv
import os
import sys
import re
import pickle
import glob
import numpy as np
import scipy.spatial.distance
from collections import defaultdict
import scipy.stats
import shutil

def intersect(a, b):
     """ return the intersection of two lists """
     ##convert all list elements to upper case
     a = [item.upper() for item in a]
     b = [item.upper() for item in b]
     return list(set(a) & set(b))

def union(a, b):
    """ return the union of two lists """
    a = [item.upper() for item in a]
    b = [item.upper() for item in b]
    return list(set(a) | set(b))

def flatten(lst):
    for elem in lst:
        if type(elem) == list:
            for item in elem:
                yield str(item)
        else:
            yield str(elem)

def calcIC(l):
    for i in l:
	rowvec = i.split('\t')
	if "A" in rowvec[0]:
	    a = rowvec[1:]
	    sumA=0
	    for e in a:
		if float(e) ==0:
		    sumA += float(e)*np.log((1e-5)/0.2953)##this is zero:/
		else:
		    sumA += float(e)*np.log(float(e)/0.2953)
	elif "C" in rowvec[0]:
	    c = rowvec[1:]
	    sumC=0
	    for e in c:
		if float(e) == 0:
		    sumC += float(e)*np.log((1e-5)/0.2047)
		else:
	            sumC += float(e)*np.log(float(e)/0.2047)
	elif "G" in rowvec[0]:
	    g = rowvec[1:]
	    sumG=0
	    for e in g:
		if float(e) == 0:
		    sumG += float(e)*np.log((1e-5)/0.2047)
		else:
		    sumG += float(e)*np.log(float(e)/0.2047)
	elif "T" in rowvec[0]:
	    t = rowvec[1:]
	    sumT=0
	    for e in t:
		if float(e)==0:
		    sumT += float(e)*np.log((1e-5)/0.2953)
		else:
		    sumT += float(e)*np.log(float(e)/0.2953)
    return sumA+sumC+sumG+sumT 

def calcED(ll1,ll2):
    ## assert matrix dimensions match                    
    if len(ll1[0]) == len(ll2[0]):
	a1 = np.array(ll1)
        a2 = np.array(ll2)
        at1 = np.transpose(a1)##nby4
        atv1 = at1.reshape(1,len(ll1[0])*4)
        at2 = np.transpose(a2)
        atv2 = at2.reshape(1,len(ll1[0])*4)
        ##check the other orientation
        af1 = np.fliplr(np.flipud(at1))
        afv1 = af1.reshape(1,len(ll1[0])*4)
        #print 'check matrix dimension', a1.shape, a2.shape
        #print atv1, atv2, afv1
        Y1 = scipy.spatial.distance.cdist(atv1,atv2,'euclidean')
        Y2 = scipy.spatial.distance.cdist(afv1,atv2,'euclidean')    
	return (Y1,Y2)
    else:
	print 'Warning: unmatched dimension'

def calcPCC(ll1, ll2):
    ## assert matrix dimensions match
    if len(ll1[0]) == len(ll2[0]):
        a1 = np.array(ll1,dtype='float64')
        a2 = np.array(ll2,dtype='float64')
        at1 = np.transpose(a1)#nby4 matrix
        atv1 = at1.reshape(len(ll1[0])*4,)
        at2 = np.transpose(a2)
        atv2 = at2.reshape(len(ll1[0])*4,)
        ##check the other orientation
        af1 = np.fliplr(np.flipud(at1))
        afv1 = af1.reshape(len(ll1[0])*4,)
        #print 'check matrix dimension', a1.shape, a2.shape
        #print at1, at2, af1
	#print atv1[0]
        #Y1 = scipy.spatial.distance.cdist(atv1,atv2,'correlation')
        #Y2 = scipy.spatial.distance.cdist(afv1,atv2,'correlation')
	Y1 = scipy.stats.pearsonr(atv1,atv2)[0]
        Y2 = scipy.stats.pearsonr(afv1,atv2)[0]
        return (Y1,Y2)
    else:
        print 'Warning: unmatched dimension'

def main(argv):
    if len(argv) < 2:
        sys.stderr.write("Usage: %s Path-to-Motif-Pool\n" % argv[0])
        return 1
    if not os.path.exists(argv[1]):
        sys.stderr.write('Error: Path %r was not found!\n' % argv[1])
	return 1
    #if not os.path.exists(argv[2]):
	#sys.stderr.write('Error: Path-2 %r was not found!\n' % argv[2])
    	#return 1

    path1 = sys.argv[1]
#    count = 0
    createSimDict = False
    if createSimDict:
        eddict = defaultdict(tuple)
        mlist = []
        ##input motif files from ENCODE
        for infile1 in glob.glob(os.path.join(path1, "*_uniprobe")):
            (PATH1,FILENAME1) = os.path.split(infile1)
	    F1 = FILENAME1.split('_uniprobe')[0]
	    try:
	        f1 = open(infile1, 'rt')
	    except IOError:
	        continue
	    str1 = f1.read()
	    list1 = str1.split('\n')
	    for infile2 in glob.glob(os.path.join(path1, "*_uniprobe")):
                (PATH2, FILENAME2) = os.path.split(infile2)
	        F2 = FILENAME2.split('_uniprobe')[0]
                print 'Comparing: '+ F1 + ',' + F2
	        #print 'list1', list1
	    	try:
	            f2 = open(infile2,'rt')
	    	except IOError:
		    continue
            	str2 = f2.read()
                list2 = str2.split('\n')
	        #print 'list2', list2
	        l1 = list1[1:]
	        l2 = list2[1:]
	        m1 = len(l1[0].split('\t'))-1
	        m2 = len(l2[0].split('\t'))-1
	        ic1 = calcIC(l1)
	        ic2 = calcIC(l2)
	        print "information content",ic1, ic2
	        offset = abs(m1-m2)
	        print offset
	        ##create a list to hold all matrix distance values 
	        dlist = []
	        ## try different offsets and both orientations
	        if offset > 0:
                    for i in xrange(offset):
		        ll1 = []
		        ll2 = []
                        if m1 > m2:
		            for j in xrange(4):
			        ll1.append(l1[j].split('\t')[1+i:1+i+min(m1,m2)])
                                ll2.append(l2[j].split('\t')[1:])
			    #l1_offset = l1i:i+min(m1,m2)]
                            #print 'l1 offset',ll1
		            #print 'l2', ll2
                        else:
		            for j in xrange(4):
			        ll1.append(l1[j].split('\t')[1:])
                                ll2.append(l2[j].split('\t')[1+i:1+i+min(m1,m2)])
           	            #print 'l1', ll1
	                    #print 'l2 offset', ll2
		        Y1,Y2 = calcPCC(ll1,ll2)
	 	        dlist.append(float(Y1))
		        dlist.append(float(Y2))
	        else:
		    ll1 = []
		    ll2 = []
		    for j in xrange(4):
		        ll1.append(l1[j].split('\t')[1:])
		        ll2.append(l2[j].split('\t')[1:])
	            Y1,Y2 = calcPCC(ll1,ll2)	 
	            dlist.append(float(Y1))
	            dlist.append(float(Y2))
	            #s = np.linalg.norm(Y)/float(a1.shape[0]+a2.shape[0])
	            #print s, list1[0], list2[0]
	            #if list1[1:] == list2[1:]:
	    	        #print 'identical' + FILENAME1 + ',' + FILENAME2
		    #print list1[0]
		    #count +=1
	    #print 'dlist', dlist 
	    ##removing identical motifs
	    	cleanMotifs = False
	    	if cleanMotifs:
	            if max(dlist) == 1.0:
		    	if not FILENAME1 == FILENAME2:
		            if '_0.90_' in FILENAME1:
			    	if not '_0.90_' in FILENAME2:
		                    print 'removing', FILENAME1
		                    os.remove(infile1)
		            else:
			    	if '_0.90_' in FILENAME2:
		                    print 'removing', FILENAME2
		                    os.remove(infile2)

	        eddict[(F1,F2)]=max(dlist)#,len(ll1[0]))
	        if not (F1,m1,ic1) in mlist:
		    mlist.append((F1,m1,ic1))
		    if not (F2,m2,ic2) in mlist:
		        mlist.append((F2,m2,ic2))

    saveDistDict = True
    if saveDistDict:
        with open('/home/xc406/code/motifPCC.pickle','w') as f:
            pickle.dump(eddict,f)
        with open("/home/xc406/code/mllist.pickle",'w') as f:
	    pickle.dump(mlist,f)

    readDistDict = True
    if readDistDict:
        with open('/home/xc406/code/motifPCCmm9.pickle','r') as f:
	    eddict = pickle.load(f)
        with open("/home/xc406/code/mllistmm9.pickle",'r') as f:
	    mllist = pickle.load(f)

    createFlist = True
    if createFlist:
        flist = []
        #mllist = mlist
        for (f1, m1, ic1) in mllist:
	    if "_0.90" in f1:
	        mname1 = f1.split("_")[:-2]
	    else:
	        mname1 = f1.split("_")
	        del mname1[1]
	    f = f1## name of motif to keep
	    m = m1## length of motif for comparison
	    ic = ic1## information content of motif for comparison
	    for (f2, m2, ic2) in mllist:
	        if "_0.90" in f2:
	            mname2 = f2.split("_")[:-2]
	        else:
	            mname2 = f2.split("_")
		    del mname2[1]
		    #if len(mname2) != 3:
		        #print "encode name", mname2
	        if (f1 != f2) and (len(intersect(mname1,mname2)) > 0):
		    print intersect(mname1,mname2)
		    if eddict[(f1,f2)] > 0.74:
		        if m2 > m:##compare motif length
		            m = m2
		            f = f2
			    ic = ic2
		        elif m2 == m:
			    if ic2 > ic:##compare information content if motif lengths equal
			        m = m2
			        f = f2
			        ic = ic2
		        print "found similar:", f1, f2, m1, m2, ic1, ic2
		        print "keep:", f,m,ic
	    if not f in flist:
	        flist.append(f)

        print len(flist)
        with open("/home/xc406/code/flist74mm9.pickle", "w") as f:
	    pickle.dump(flist, f)
        #with open("/home/xc406/code/flist.pickle", "r") as f:
	    #flist = pickle.load(f)

    moveFiles = True
    if moveFiles:
        for f in flist:
	    f += "_uniprobe"
	    #print path1+'keep'
	    shutil.copy(os.path.join(path1,f), path1+'keep74')	

    writeMotifMat = True
###write distance matrix output
    if writeMotifMat:
        ofile = open('/scratch/xc406/hg19fimo/MotifPccMatmm9all', 'w')
        writer = csv.writer(ofile,delimiter = '\t')

        mlist = []
        for m1,m2 in eddict:
            if not m1 in mlist:
                mlist.append(m1)
                if not m2 in mlist:
                    mlist.append(m2)
        mlist.sort()
        flist_name = []
        for item in mlist:
	    ## exclude unannotated motifs
            if not "Noname" in item:
                flist_name.append(item)
        flist_name.sort()
        writer.writerows([list(flatten(['motif', flist_name]))])

        flist = mlist

        for i in xrange(len(flist)):
            if not "Noname" in flist[i]:
                line = []
                line.append(flist[i])
                for j in xrange(len(flist)):
                    if not "Noname" in flist[j]:
                        try:
                            line.append(eddict[(flist[i],flist[j])])
                        except KeyError:
                            line.append('NA')
                            print 'missing pairwise comparison', flist[i], mlist[j]
                if len(line) == len(flist_name)+1:
                    writer.writerows([line])
                else:
                    print 'writer error!'

        ofile.close()
	
if __name__=='__main__':
    sys.exit(main(sys.argv))

