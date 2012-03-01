"""
released under the terms of the LGPL
copyright ross lazarus August 2007 
for the rgenetics project

Special galaxy tool for the camp2007 data
Allows grabbing arbitrary columns from an arbitrary region

Needs a mongo results file in the location hardwired below or could be passed in as
a library parameter - but this file must have a very specific structure
rs chrom offset float1...floatn

called as
    <command interpreter="python">
        camp550KV3.py $range $name $out_file1
    </command>

cols is a delimited list of chosen column names for the subset
r is a ucsc location region pasted into the tool

"""


import sys

DATA_FILE = '/usr/local/galaxy/data/camp550KV3/550KV3_PBAT_auto.ped'
MAP_FILE = '/usr/local/galaxy/data/camp550KV3/camp2007.map'

### Sanity check the arguments
if len(sys.argv) != 3: 
    print >> sys.stdout, '##!camp550KV3: Expected 2 params in sys.argv, got %d (%s)' % (len(sys.argv), sys.argv)
    sys.exit(1)

### Figure out what genomic region we are interested in
region = sys.argv[1]
try:
    c,rest = region.split(':')
    c = int(c.replace('chr', ''))   
    rest = rest.replace(',', '')
    spos, epos = rest.split('-')
    spos = int(spos)
    epos = int(epos)
except:
    print >> sys.stdout, '##!camp550KV3: Unable to parse region "%s" - MUST look like "chr8:10,000-100,000' % (region)
    sys.exit(1)

### Figure out which markers are in this region
mfile = open(MAP_FILE, 'r')
markers = []
for line in mfile:
    line = line.strip()
    if not line: continue
    chrom, snp, _, abspos = line.split()
    abspos = int(abspos)
    chrom = int(chrom)
    if chrom == c and spos <= abspos <= epos:
        markers.append(snp)
mfile.close()

### Now parse the first line of the ped file to determine which columns
### have the markers we want ...
try:
    dfile = open(DATA_FILE, 'r')
except: # bad input file name?
    print >> sys.stdout, '##!camp550KV3 unable to open file %s' % (DATA_FILE)
    raise
    sys.exit(1)
header = dfile.readline().strip()
snps = header.split()
snpcols = []
for index, snp in enumerate(snps):
    if snp in markers:
        snpcols.append(index)

### ... and then parse the rest of the ped file to pull out
### the genotypes for all subjects for those markers
res = []
for line in dfile:
    line = line.strip()
    if not line: continue
    fields = line.split()
    preamble = fields[:6]
    genotypes = []
    for snpcol in snpcols:
        genotypes.append('%s\t%s' % (fields[snpcol+6], fields[snpcol+7]))
    res.append('%s\t%s' % ('\t'.join(preamble), '\t'.join(genotypes)))
dfile.close()    

### Write the output header ...
ofname = sys.argv[2]
ofile = file(ofname, 'w')
ofile.write('%s\n' % '\t'.join(markers))

### ... And the data
for line in res:
    ofile.write('%s\n' % (line))
ofile.close()


