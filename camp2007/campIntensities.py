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
        campRS.py /usr/local/galaxy/data/camp2007/camp2007.xls $cols "$r" $name $out_file1 "$rslist"
    </command>

cols is a delimited list of chosen column names for the subset
r is a ucsc location region pasted into the tool

"""


import sys          

if len(sys.argv) < 8: 	
  print >> sys.stdout, '##!expected 8 params in sys.argv, got %d - %s' % (len(sys.argv),sys.argv)
  sys.exit(1)
# quick and dirty for galaxy - we always get something for each parameter
fname = sys.argv[1]
wewant = sys.argv[2].split(',')
region = sys.argv[3].strip()
name = sys.argv[4]
ofname = sys.argv[5]
listrs = sys.argv[6]
map_path = sys.argv[7]
if len(wewant) == 0: # no columns selected?
  print >> sys.stdout, '##!campIntensities - no columns selected - cannot run'
  sys.exit(1)
try:
  f = open(fname,'r',2**30)
except: # bad input file name?
  print >> sys.stdout, '##!campIntensities unable to open file %s' % (fname)
  sys.exit(1)
rslist = []
rsdict = {}
c = None
spos = 0
epos = 0
if region > '':
    try: # TODO make a regexp?
        c,rest = region.split(':')
        c = c.replace('chr','')
        try:
            c = int(c)
        except:
            pass
        rest = rest.replace(',','') # remove commas
        spos,epos = rest.split('-')
        spos = int(spos)
        epos = int(epos)
        print >> sys.stdout, '##campIntensities parsing chrom %s from %d to %d' % (c,spos,epos)
    except:
        print >> sys.stdout, '##!campIntensities unable to parse region %s - MUST look like "chr8:10,000-100,000' % (region)
        sys.exit(1)
    mfile = open(map_path, 'r')
    for line in mfile:
        if not line:
            continue
        chrom, snp, gpos, abspos = line.strip().split()
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
    rslist = rslist.replace('__cn__',' ').split() # galaxy replaces newlines with __cn__ - go figure
rsdict = dict(zip(rslist,rslist)) # for quick lookups
print >> sys.stdout, '##campIntensities looking for %d rs (%s)' % (len(rslist),rslist[:5])
res = []
head = "              "
while head[:9] <> 'Sample ID':
      head = f.readline()
cnames = head.strip().split('\t')
wewant = map(int,wewant)
minline = max(wewant) # largest column 
for n,l in enumerate(f):
  ll = l.strip().split()
  if len(ll) >= minline:
    thisrs = ll[1]
    if rsdict.get(thisrs,None):
      row = [ll[x] for x in wewant]
      res.append(row) # subset of columns plus sort keys
o = file(ofname,'w')
res.sort()
res = ['%s\n' % '\t'.join(x) for x in res] # turn into tab delim string without the prefix chrom offset
print >> sys.stdout, '##campIntensities selected and returning %d data rows' % (len(res))
head = [cnames[x] for x in wewant] # ah, list comprehensions - list of needed column names
o.write('%s\n' % '\t'.join(head)) # header row for output
o.write(''.join(res))
o.close()
f.close()    


