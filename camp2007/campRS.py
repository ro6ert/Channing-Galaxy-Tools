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

if len(sys.argv) < 7: 
  print >> sys.stdout, '##!expected 7 params in sys.argv, got %d - %s' % (len(sys.argv),sys.argv)
  sys.exit(1)
# quick and dirty for galaxy - we always get something for each parameter
fname = sys.argv[1]
wewant = sys.argv[2].split(',')
region = sys.argv[3].strip()
name = sys.argv[4]
ofname = sys.argv[5]
listrs = sys.argv[6]
if len(wewant) == 0: # no columns selected?
  print >> sys.stdout, '##!campRS - no columns selected - cannot run'
  sys.exit(1)
try:
  f = open(fname,'r')
except: # bad input file name?
  print >> sys.stdout, '##!campRS unable to open file %s' % (fname)
  sys.exit(1)
rslist = []
rsdict = {}
c = None
galsep1 = '__cr__'
galsep2 = '__cn__' # updated feb 2010
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
        print >> sys.stdout, '##campRS parsing chrom %s from %d to %d' % (c,spos,epos)
    except:
        print >> sys.stdout, '##!campRS unable to parse region %s - MUST look like "chr8:10,000-100,000' % (region)
        sys.exit(1)
else:
    listrs = listrs.replace(galsep1,' ') # in case only __cr__ one space or two doesn't matter
    rslist = listrs.replace(galsep2,' ').lower().split() # galaxy replaces newlines with XX - go figure
    rsdict = dict(zip(rslist,rslist)) # for quick lookups
    print >> sys.stdout, '##campRS looking for %d rs (%s)' % (len(rslist),rslist[:5])
res = []
cnames = f.next().strip().split() # column titles for output
wewant = map(int,wewant)

for n,l in enumerate(f):
  ll = l.strip().split()
  thisrs = ll[0]
  try:
      thisc = int(ll[1])
  except:
      thisc = ll[1]
  try:
      thispos = int(ll[2])
  except:
      print >> sys.stdout, '##campRS dud offset %s at %d' % (ll[:5],n)      
      thispos = 0
  if rsdict.get(thisrs,None) or ((thisc == c) and (thispos >= spos) and (thispos <= epos))  :
     row = [ll[x] for x in wewant]
     row.insert(0,thispos)
     row.insert(0,thisc)
     res.append(row) # subset of columns plus sort keys
o = file(ofname,'w')
res.sort()
res = ['%s\n' % '\t'.join(x[2:]) for x in res] # turn into tab delim string without the prefix chrom offset
print >> sys.stdout, '##campRS selected and returning %d data rows' % (len(res))
head = [cnames[x] for x in wewant] # ah, list comprehensions - list of needed column names
o.write('%s\n' % '\t'.join(head)) # header row for output
o.write(''.join(res))
o.close()
f.close()    


