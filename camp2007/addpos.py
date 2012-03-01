# quick hack to add coordinates to all the pbat
# run results from amy for the CAMP 2007 illumina 
# run - sorts and inserts chrom offset at start of line
# actually after rs 
# this file is used by the special camp subsetting tool
# that can subset by genomic region and columns 
# for analyses by the fellows
#
# ross lazarus
# 
# map file
m = open('camp2007.map','r')
map = [x.strip().split() for x in m]
rsdict = {}
for x in map:
  c,rs,gpos,ppos = x
  rsdict[rs] = (c,ppos)
# all the results file 
i = open('Aug2007trackfile.txt','r')
head = i.next().strip().split()
head.insert(1,'Offset')
head.insert(1,'Chrom')
# add our new cols
o = open('/usr/local/galaxy/data/camp2007/camp2007.xls','w')
# ugh hardwired too
s = '\t'.join(head)
o.write('%s\n' % s) # write the new header
res = [] # sort!
for n,l in enumerate(i):
   ll = l.strip().split()
   rs = ll[0]
   c,ppos = rsdict.get(rs)
   ll.insert(1,int(ppos))
   ll.insert(1,int(c))
   res.append(ll)
   if (n+1) % 100000 == 0:
      print 'at line',n,ll
print 'sorting'
res.sort()
print 'fixing'
for n in range(len(res)):
   r = res[n]
   r[2] = '%d' % r[2]
   r[1] = '%d' % r[1]
   o.write('%s\n' % '\t'.join(r))
   if (n+1) % 100000 == 0:
        print 'at line',n,r
o.close()
i.close()
m.close()

