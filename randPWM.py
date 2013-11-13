import numpy as np
from Bio import SeqIO
from Bio import motifs

handle = open("/home/xc406/data/hg19randomnorepeat.fa", "rU")
records = list(SeqIO.parse(handle, "fasta"))

instances = []
for record in records:
    instances.append(record.seq)
m = motifs.create(instances)
pwm = m.counts.normalize(pseudocounts={'A':0.2953, 'C': 0.2047, 'G': 0.2047, 'T': 0.2953})

na ='\t'.join(list(map(str,pwm['A'])))
nc ='\t'.join(list(map(str,pwm['C'])))
ng ='\t'.join(list(map(str,pwm['G'])))
nt ='\t'.join(list(map(str,pwm['T'])))
r1 = '\t'.join(['A:',na])
r2 = '\t'.join(['C:',nc])
r3 = '\t'.join(['G:',ng])
r4 = '\t'.join(['T:',nt])
r = '\n'.join(['M0000_0.80',r1,r2,r3,r4])

with open('randompwm.txt','w') as f:
    f.write(r)

with open('hg19chromdict.txt', 'r') as f:
        chromdict = pickle.load(f)

chromdict.rand
r = np.random.rand(10,4)
new_matrix = np.zeros((10,4))
row_sums = r.sum(axis=1)
for i, (row,row_sum) in enumerate(zip(r,row_sums)):
    new_matrix[i,:] = row/row_sum

##transpose
