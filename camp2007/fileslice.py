#!/usr/local/bin/python2.4
# quick and dirty file slicer for vince
# copyright ross lazarus 2008
# released under the LGPL
# for the rgenetics project
# added headerskiplines param for Jen-wah

import sys,os,time

def fslice(fname='',n=5614660,prog='?',hsl=0):
   """make sub files from fname of n lines each
   skip hsl lines 
   """
   f = file(fname,'r')
   newfname = os.path.split(fname)[-1]
   outf = None
   started = time.time()
   for i,row in enumerate(f):
       if i >= hsl: # do nothing until then
           if ((i-hsl) % n == 0): # time to write if not first time
               if outf: # not first time
                 outf.close()
                 outf = None
           if outf == None:
              newf = '%d_%s' % (i-hsl,newfname)
              outf = file(newf,'w')
              dur = time.time()-started
              if i > hsl:
                 print '%s: %s opened, duration = %e secs, %e rows/sec' % (prog, newf, dur, i/dur)
           outf.write(row)
        else:
            if i % 1000000 == 0:
                print 'at line %d of %d lines to be skipped' % (i,hsl)
   if outf:
       outf.close()

if __name__ == "__main__":
    progname = os.path.split(sys.argv[0])[-1]
    if len(sys.argv) > 3:
        n = int(sys.argv[1])
        fname = sys.argv[2]
        hsl = int(sys.argv[3])
        fslice(fname,n,progname,hsl)
    else:
        print "%s needs 3 parameters" % progname
        print "the number of lines per slice, the input file name, and header rows to skip (0=none)" 
        sys.exit(1)

