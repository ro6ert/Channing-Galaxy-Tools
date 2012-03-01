# for rgenetics - lped to fbat
# recode to numeric fbat version
# much slower so best to always
# use numeric alleles internally

import sys,os

myname = os.path.split(sys.argv[0])[-1]


def fbater(pedname=''):
  """convert linkage ped/map to fbat"""  
  recode={'A':'1','C':'2','G':'3','T':'4','N':'0','0':'0','1':'1','2':'2','3':'3','4':'4'}  
  inmap = '%s.map' % pedname
  try:
    mf = file(inmap,'r')
  except:
    sys.stderr.write('%s cannot open inmap file %s - do you have permission?\n' % (myname,inmap))
    sys.exit(1)
  try:
    rsl = [x.split()[1] for x in mf]
  except:
      sys.stderr.write('## cannot parse %s' % inmap)
      sys.exit(1)   
  head = ' '.join(rsl) # list of rs numbers
  # TODO add anno to rs but fbat will prolly barf?
  outped = '%sfbat.ped' % pedname
  inped = '%s.ped' % pedname
  pedf = file(inped,'r')
  o = file(outped,'w',2**30)
  o.write(head)
  o.write('\n')
  for i,row in enumerate(pedf):
    if (i+1) % 1000 == 0:
        print '%s:#%d' % (inped,i+1)
    if i == 0:
       lrow = row.split()
       try:
          x = [int(x) for x in lrow[10:20]] # look for non numeric codes
       except:
          dorecode = 1
       print 'dorecode=',dorecode,'lrow=',lrow[:80]
    if dorecode:
        lrow = row.split()
        p = lrow[:6]
        g = lrow[6:]
        gc = [recode.get(x,'0') for x in g]
        lrow = p+gc
        row = '%s\n' % ' '.join(lrow)
        if i < 5:
             print '%d: lrow=%s' % (i,lrow[:50])
    o.write(row)
  o.close()


def main():
  """call fbater
  need to work with rgenetics composite datatypes
  so in and out are html files with data in extrafiles path
  """
  nparm = 3
  if len(sys.argv) < nparm:
    sys.stderr.write('## %s called with %s - needs %d parameters \n' % (myname,sys.argv,nparm))
    sys.exit(1)
  inpedname = sys.argv[1]
  outname = sys.argv[2]
  fbater(inpedname,inoutname)


if __name__ == "__main__":
   main()
