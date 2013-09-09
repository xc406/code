import sys
import os
import itertools
import csv
import pickle
#import operator

def main(argv):
    if len(argv) < 2:
        sys.stderr.write("Usage: %s gff_file\n" % argv[0])
        return 1
    if not os.path.isfile(argv[1]):
        sys.stderr.write('Error: gff_file %r was not found!\n' % argv[1])
        return 1
    
    infile = sys.argv[1]

    (path,fname) = os.path.split(infile)

    ifile = open(infile,'rt')
    reader = csv.reader(ifile, delimiter = '\t')

    ofile = open(os.path.join(path, fname + '_clean'), 'w')
    writer = csv.writer(ofile, delimiter = '\t')
	
#sortedlist = sorted(reader, key=operator.itemgetter(3), reverse=False)	    

    tgtlist = list([(((row[-1].split(';')[0]).split(' ')[1]).split('_')[0],int(row[3]),int(row[4])-int(row[3])+1,float(row[5]),row[6],row[1],row[0]) for row in reader])##list of (target,start,mlen,pval,strand,motif,chr)
    tgtdict = {}
    for t,start,mlen,p,strand,m,chr in tgtlist:
	if not t in tgtdict:
            tgtdict[t]= [(m,chr,start,mlen,p,strand)]
        else:
            tgtdict[t].append((m,chr,start,mlen,p,strand))
    ##store all the motif hits for each target gene in a dictionary
    #print tgtdict
    targetlist = []
    hitslist = []
    newtgtdict = {}
    ##remove all identicals and overlaps
    while tgtdict:
	target,hits = tgtdict.popitem()
	#hits_comp = list(hits)
	hits_new = list(hits)
	for h1 in hits:
	    if h1 in hits_new: 
	        hits_comp=list(hits_new)
	        hits_comp.remove(h1)
	        done = False##flag to break out of the h2 loop if h1 is removed
	        for h2 in hits_comp:
		    if done:
		        break
		    if h2 in hits_new:
		        if abs(h2[2]-h1[2]) < min(h2[3],h1[3])/2:##compare motif hits of the same target when they have over half of the sites overlap
		            if h2[4] < h1[4]:##compare p vals and keep the one with smaller p val
			        hits_new.remove(h1)
			        done = True
		            else:
			        hits_new.remove(h2) 
	
	#targetlist.append(target)
	#hitslist.append(hits_new)
	newtgtdict[target] = hits_new
	#print hitslist	
    #print len(targetlist), len(hitslist)
    #newtgtdict=dict(itertools.izip(targetlist,hitslist))
    #print newtgtdict		      
    #with open('tgtdict.txt','w') as f:
        #pickle.dump(newtgtdict,f)	
        
    while newtgtdict:
	target,hits = newtgtdict.popitem()
	for h in hits:
	    ##gff_outfile[chr,motif_ID,'motif',start,end,pval,strand,.,group(gene_id:gene_name_chr_start)]
	    newrow = list([h[1],h[0],'motif',h[2],h[2]+h[3]-1,h[4],h[5],'.','gene_id '+ target + '_' + h[1] + '_' + str(h[2]) + '; seq ']) 
	    writer.writerows([newrow])
    ofile.close()
	

if __name__=='__main__':
    sys.exit(main(sys.argv))



