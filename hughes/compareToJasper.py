import csv

ifile = open('/Users/xclair/Documents/google-python-exercises/hg','r')
reader = csv.reader(ifile, delimiter = '\t')
ofile = csv.writer(open('/Users/xclair/Documents/google-python-exercises/NM_001204184','w'), delimiter = '\t')

for row in reader:
    if row[2] == 'NM_001204184':
##    if 'TFAP2A' in row[1]:
        print row
        ofile.writerows([row])

