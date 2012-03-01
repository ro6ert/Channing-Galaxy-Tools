# join 2 eigenstrat geno files
# expect path to two sets to join (must have .map, .ind and .eigenstratgeno)
# and a path for the new set
# maps are checked - cannot join if not same
# ind files are concatenated
# eigenstratgeno files are joined row-wise

import os, sys, shutil

def mapsMatch(m1,m2):
   """
   maps must match for joining to work
   """
   same = True
   f1 = file(m1,'r').readlines()
   f2 = file(m2,'r').readlines()
   for i, row in enumerate(f1):
       row = row.strip().split()
       row2 = f2[i].strip().split()
       if row[0] <> row2[0]:
	    same = False
            break
   return same 


def joinRows(r1,r2,outfname):
   """
   geno codes must be appended row (marker) wise in subject (.ind file) order
   """
   outf = open(outfname,'w')
   f1 = file(r1,'r')
   f2 = file(r2,'r')
   for row1 in f1:
       if row1.strip() > '':
          row2 = f2.next()
          outf.write('%s%s\n' % (row1.strip(),row2.strip()))
   outf.close()


def joinInds(r1,r2,outfname):   
   """
   individual data must be appended in same order as genos are being added
   """
   outf = open(outfname,'w')
   f1 = file(r1,'r')
   f2 = file(r2,'r')
   for row1 in f1:
       outf.write('%s\n' % (row1.strip()))
   for row1 in f2:
       outf.write('%s\n' % (row1.strip()))
   outf.close()


def dojoin(ipath1,ipath2,opath):
    """ provide both input file names up to extension
    and outfile path including name up to extension
    """
    r1 = '%s.map' % ipath1
    r2 = '%s.map' % ipath2
    if not mapsMatch(r1,r2):
         print '### maps %s and %s do not match' % (r1,r2)
         sys.exit(1)
    outpath = '%s.map' % opath
    shutil.copyfile(r1,outpath)
    r1 = '%s.eigenstratgeno' % ipath1
    r2 = '%s.eigenstratgeno' % ipath2
    outpath = '%s.eigenstratgeno' % opath
    joinRows(r1,r2,outpath)
    outpath = '%s.ind' % opath
    r1 = '%s.ind' % ipath1
    r2 = '%s.ind' % ipath2
    joinInds(r1,r2,outpath)

if __name__ == "__main__":
    progname = os.path.split(sys.argv[0])[-1]
    if len(sys.argv) >= 4:
	dojoin(sys.argv[1],sys.argv[2],sys.argv[3])
    else:
        print 'provide all but the extension for 2 sets of eigenstratgeno files and for the output'
        print """eg python %s python rgJoinEig.py 
/opt/galaxy/database/files/001/dataset_1023_files/Pollak_infiniumblack_remocontrnonwh 
/opt/galaxy/database/files/001/dataset_1042_files/Pollak_infiniumwhites_remcontrol  
/opt/galaxy/tool-data/rg/library/eigenstratgeno/Pollack_infinium_allCases""" % progname
        sys.exit(1)
