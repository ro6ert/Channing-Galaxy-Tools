# split a bed file into chr files
# should preserve order but no sort as at November 2010
# Ross Lazarus
# Copyright 2010. All rights reserved.
# Distributed under the terms of the LGPL.

import sys,os

if len(sys.argv) == 1:
   print 'Need an input bed file name to split into chromosomes'
   sys.exit(1)
fname = sys.argv[1]
assert os.path.isfile(fname), '## Input error: "%s" is not a bed filename we seem to be able to split into chromosome files' % fname
fouts = {}
seen = {}
f = file(fname,'r')
froot = os.path.splitext(fname)[0]
for i,row in enumerate(f):
    c = row.split('\t')[0]
    seen.setdefault(c,0)
    seen[c] += 1
    if seen[c] == 1: # first time
        print 'seeing %s' % c
        cfname = '%s_%s.bed' % (froot,c)
	fouts[c] = file(cfname,'w') 
    fouts[c].write(row) # should preserve order seen
k = fouts.keys()
k.sort()
for c in k:
    fouts[c].close()
    print '# wrote %d rows to %s' % (seen[c],c)


