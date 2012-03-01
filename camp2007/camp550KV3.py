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
        camp550KV3.py '$region' '$rslist' '$title' $out_file1
    </command>

cols is a delimited list of chosen column names for the subset
r is a ucsc location region pasted into the tool

"""


import sys,os

###DATA_FILE = '/usr/local/galaxy/data/camp550KV3/550KV3_PBAT_auto.ped'
###MAP_FILE = '/usr/local/galaxy/data/camp550KV3/camp2007.map'
DATA_FILE = '/share/shared/data/camp550KV3/550KV3_PBAT_auto.ped'
MAP_FILE = '/share/shared/data/camp550KV3/camp2007.map'

progname = os.path.split(sys.argv[0])[1]

### Sanity check the arguments
if len(sys.argv) != 5: 
    print >> sys.stdout, '''##! Rgenetics Genomic Region Subsetter:
    Expected 5 params in sys.argv, got %d (%s)''' % (len(sys.argv), sys.argv)
    sys.exit(1)

### Figure out what genomic region we are interested in
region = sys.argv[1].strip()
listrs = sys.argv[2].strip()
title = sys.argv[3].replace(' ','_')
ofname = sys.argv[4]
c = None
spos = epos = 0
rsdict = {}
rslist = []
if region > '':
    try:
        c,rest = region.split(':')
        c = c.replace('chr', '')
        try:
            c = int(c)
        except:
            pass
        rest = rest.replace(',', '')
        spos, epos = rest.split('-')
        spos = int(spos)
        epos = int(epos)
    except:
        print >> sys.stdout, '''##! Rgenetics Genomic Region Subsetter:
        Unable to parse region "%s" - MUST look like "chr8:10,000-100,000''' % (region)
        sys.exit(1)
    ### Figure out which markers are in this region
    mfile = open(MAP_FILE, 'r')
    for line in mfile:
        line = line.strip()
        if not line:
            continue
        chrom, snp, _, abspos = line.split()
        abspos = int(abspos)
        try:
            chrom = int(chrom)
        except:
            pass
        if (chrom == c and spos <= abspos <= epos):
            rslist.append(snp)
    mfile.close()
else:
    rslist = listrs.lower().replace('__cr__',' ')
    rslist = rslist.replace('__cn__',' ').split() # galaxy replaces newlines with XX - go figure
    
rsdict = dict(zip(rslist,rslist)) # for quick lookups
print >> sys.stdout, '## Rgenetics %s looking for %d rs (%s)' % (progname, len(rslist),rslist[:5])


### Now parse the first line of the file to determine which columns
### have the snps we want ...
try:
    dfile = open(DATA_FILE, 'r')
except: # bad input file name?
    print >> sys.stdout, '##! Rgenetics %s unable to open file %s' % (progname,DATA_FILE)
    raise
    sys.exit(1)
header = dfile.readline().strip()
snps = header.split()
snpcols = []
omlist = []
for index, snp in enumerate(snps):
    if rsdict.get(snp,None):
        snpcols.append(6+(2*index)) # ignore first 6 pedigree fields
        omlist.append(snp) # make sure we get the right output order
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
        genotypes.append('%s %s' % (fields[snpcol], fields[snpcol+1])) # pair starting here
    res.append('%s %s' % (' '.join(preamble), ' '.join(genotypes)))
dfile.close()    

### Write the output header ...

ofile = file(ofname, 'w')
ofile.write('%s\n' % ' '.join(omlist))

### ... And the data
for line in res:
    ofile.write('%s\n' % (line))
ofile.close()


