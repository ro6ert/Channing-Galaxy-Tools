"""
look for bad q/s lengths - failed grooming on ENCODE fastq
from ucsc. Ann and friends will investigate

Looks like the whole batch are affected.
Ugly solution - assume it's the only error and fix it

Valid fastq format data look like this.

Mon,Nov 15 at 7:59am head *.fastq
@HWI-EAS68_6_FC206E3_1_1_118_667                                                                                                                                                   
GAAATTATTTTTTCCGAATTGAAGATGAAAATA
+                                                                                                                                                                                  
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIA
@HWI-EAS68_6_FC206E3_1_1_190_636                                                                                                                                                   
GAGAACACATTTTCTCACTGTTGAGCTAATAAT
+                                                                                                                                                                                  
IIIIIIIIIIIIIIIIIIIIIIIIIIIIII1CI
@HWI-EAS68_6_FC206E3_1_1_180_558                                                                                                                                                   
GTTTCCTAAATTGTAAGTTGAAAGATTTAGAAG

Written to test some ENCODE fastq files downloaded from the UCSC repo

Copyright Ross Lazarus at gmail dot com 2010
Source code licensed under the terms of the LGPL
enjoy...

"""
import os,sys,gzip,optparse

verbose = True

class fastQ:
    
    
    def __init__(self,argv,opts):
        """ simple representation of 4 rows of fastq data
        """
        self.recsPerRec = opts.recPerRec
        we = argv[0]
        self.verbose = opts.verbose
        self.ncat1 = 1 # seq row
        self.ncat2 = 3 # qual row
        self.progname = os.path.split(we)[-1] or 'fastqClassy'
        self.fname = opts.input
        self.fixifThing = opts.fixif # single characters to ignore if only 1 byte seq/qual length difference 
        self.outfname = opts.output # fixed version - trailing bogus (eg !) character removed if lengths same
        if verbose:
            print 'hello from %s with fname=%s and fixifThing=%s' % (self.progname,self.fname,self.fixifThing)
        if self.outfname:
            if opts.gzipout:
                self.outf = gzip.open(self.outfname,'wb')
            else:
                self.outf = open(self.outfname,'w')
        self.cats=opts.cats or ('ignore','seq','ignore','qual')
        # default is fastq
        self.logf = open(opts.logfname,'w')
        ext = os.path.splitext(self.fname)[-1]
        if ext == '.gz':
            self.dataf = gzip.open(self.fname,'r')
        else:
            self.dataf = open(self.fname,'r')
        
    def shutdown(self):
        """
        """
        if self.outfname:
            self.outf.close()    
        self.logf.close()
        self.dataf.close()
    
    def processRec(self,thisrec):
        """ parse a multiline fastq sequence record
        """
        for ncat,row in enumerate(thisrec):        
            if ncat == self.ncat1:
                slen = rlen
                srow = row
                qlen = None # coming soon
            elif ncat == self.ncat2:
                qlen = rlen
                qrow = row
                if slen <> qlen:
                    if abs(slen-qlen) == 1: # exactly one different check encode kludge repair case
                        if slen > qlen:
                            testme = srow[-1]
                        else:
                            testme = qrow[-1]
                        if testme in self.fixifThing:
                            # allow to pass
                            if slen > qlen:
                                thisrec[self.ncat1] = thisrec[self.ncat1][:-1] # trim
                            else:
                                thisrec[self.ncat2] = thisrec[self.ncat2][:-1] # trim
                            writerec = True
                            s= 'Repaired seq %s qual: %s\n' % (srow,qrow) #
                            self.logf.write(s)
                            self.fixed += 1
                            self.duds += 1
                        else:
                            s= '### seq row %d has %d seq but row %d has %d qual: %s\n' % (i-1,slen,i+1,qlen,row) # python 0=row 1
                            self.logf.write(s)
                            self.duds += 1
                    else:        
                        s= '### seq row %d has %d seq but row %d has %d qual: %s\n' % (i-1,slen,i+1,qlen,row) # python 0=row 1
                        self.logf.write(s)
                        self.duds += 1
                else:
                    writerec = True
                if writerec and self.outfname:
                    self.outf.write('\n'.join(thisrec))
                    self.outf.write('\n')
        
    def check(self):
        """
        fastq specific - should be abstract
        
        """
        slen = None
        qlen = None
        i = 0
        duds = 0
        fixed = 0
        done = False
        thisrec = []
        for i,row in enumerate(self.dataf):
            row = row.strip()
            rlen = len(row)
            writerec = False
            ncat = (i % self.recsPerRec) # fastq has 4 rows per sequence
            if ncat == 0: # init or time to process a record
                if len(thisrec) > 0:
                    self.processRec(thisrec)
                thisrec = []
            thisrec.append(row)
            if self.verbose and i % 10000000 == 0:
                s = 'Processing row %d:fixed %d of %d dud sequences found so far' % (i,self.fixed,self.duds)
                print s
        if len(thisrec) > 0: # last one
            self.processRec(thisrec)
        if self.duds > 0 or self.verbose:
            s = '%s found %d fastq sequences in error in %s. Fixed %d\n' % (self.progname,self.duds,self.fname,self.fixed)
            self.logf.write(s)
         
def parseFastq(argv=[]):
    """
    """     
    op = optparse.OptionParser()
    op.add_option('-i', '--input', default=None)
    op.add_option('-o', '--output', default=None)
    op.add_option('-l', '--logfname', default=None)
    op.add_option('-f', '--fixif', action="append",type="string")
    op.add_option('-v', '--verbose', action="store_true", dest="verbose", default=False)
    op.add_option('-z', '--gzipout', action="store_true", dest="gzipout", default=False)
    op.add_option('-c','--cats', action="append",type="string")
    opts, args = op.parse_args()
    opts.recsPerRec = 4 # for fastq
    assert opts.input <> None
    assert os.path.isfile(opts.input)
    fname = opts.input
    if not opts.logfname:
        opts.logfname = '%s_parseFastq.log' % fname    
    if opts.cats == None:
        opts.cats = ['ignore','seq','ignore','qual'] 
    if opts.fixif == None:
        opts.fixif = ['!']     
    assert opts.input, 'Please supply a valid fastq file (.fastq.gz or .fastq) as the first parameter'
    fname = opts.input
    assert os.path.isfile(fname), 'Please supply a valid fastq file (.fastq.gz or .fastq) as the first parameter'
    fixifThing=[] # the encode data I've looked at have a ! 
                  # added to qual strings starting about 1/3 of the way in to a 7GB file.
    f = fastQ(argv=sys.argv,opts=opts)
    f.check()
    f.shutdown()
        
if __name__ == "__main__":
    parseFastq()
