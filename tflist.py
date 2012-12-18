import csv
import itertools

ifile = open('/Users/xclair/Documents/google-python-exercises/mm9','r')
ifile1 = open('/Users/xclair/Documents/google-python-exercises/mm9TF_Information.txt', 'rt')
##ifile3 = open('/Users/xclair/Documents/google-python-exercises/motif_names.txt','rt')

reader = csv.reader(ifile, delimiter = '\t')
reader1 = csv.reader(ifile1, delimiter = '\t')

mylist = list([row[1] for row in reader]) 

tfs = []
for t in mylist:
        if ',' in t:
                tfs.append(t.split(', '))
        else:
                tfs.append(t)                


tfs.sort()

def flatten(lst):
	for elem in lst:
		if type(elem) == list:
			for item in elem:
				yield str(item)
		else:
			yield str(elem)

tfflattened = list(flatten(tfs))
##tfflat = list(itertools.chain.from_iterable(tf))
tf = []
for i in tfflattened:
	if not i in tf:
		tf.append(i)
tf.sort()

ofile = open('/Users/xclair/Documents/google-python-exercises/mm9tflist.txt','wt')

for item in tf:
        ofile.write(item + '\n')
ofile.close()

print len(tf)








