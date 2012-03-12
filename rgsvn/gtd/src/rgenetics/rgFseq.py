import os,string,sys,optparse,shutil,tempfile
from subprocess import Popen
from rgutils import galhtmlprefix,galhtmlpostfix,galhtmlattr,timenow 
progname = os.path.split(sys.argv[0])[1]

"""
fix mate pair information
Copyright ross lazarus oct 31 (happy halloween!)
All rights reserved
released for rgenetics under the LGPL
"""



class fixMate():
    """
    classy!
    kind of silly to make a class here.
    """
    
    def __init__(self,opts=None,cl=[],tidy=True):
        """
        """
        self.opts = opts
        self.tidy = tidy
        self.cl = ' '.join(cl) # ready for the htmlfile output
        self.delme = []
        killme = string.punctuation + string.whitespace
        trantab = string.maketrans(killme,'_'*len(killme))
        self.title = self.opts.title.translate(trantab)
        fd,self.tlogname = tempfile.mkstemp(dir=self.opts.tmp_dir,suffix='rgFixMate.log')
        self.tlog = open(self.tlogname,'w')
    

    def run(self,cl=None,redir=True):
        assert cl <> None
        fd,templog = tempfile.mkstemp(dir=self.opts.tmp_dir,suffix='rgSortBamRun.log')
        tlf = open(templog,'w')
        if redir:
           process = Popen(' '.join(cl), shell=True, stderr=tlf, stdout=tlf, cwd=self.opts.tmp_dir)
        else:
           process = Popen(' '.join(cl), shell=True, cwd=self.opts.tmp_dir)
        rval = process.wait()
        tlf.close()
        tlogs = ''.join(open(templog,'r').readlines())
        if len(tlogs) > 1:
            s = '## executing %s returned status %d and log (stdout/stderr) records: \n%s\n' % (' '.join(cl),rval,tlogs)
        else:
            s = '## executing %s returned status %d. Nothing appeared on stderr/stdout\n' % (' '.join(cl),rval)
        os.unlink(templog) # always
        self.tlog.write(s)

        
    def bamToSam(self,infile=None):
        """
        use samtools view to convert bam to sam
        """
        fd,tempsam = tempfile.mkstemp(dir=self.opts.tmp_dir,suffix='rgSortBamTemp.sam')
        cl = ['samtools view -h -o',tempsam,infile]
        self.run(cl)
        return tempsam

    def samToBam(self,infile=None):
        """
        use samtools view to convert sam to bam
        """
        fd,tempbam = tempfile.mkstemp(dir=self.opts.tmp_dir,suffix='rgSortBamTemp.sam')
        cl = ['samtools view -h -b -S -o',tempbam,infile] 
        self.run(cl)
        return tempbam
    
    def doFix(self):
        """
        """
        if self.opts.newformat == 'sam':
            fd,fixedbam = tempfile.mkstemp(dir=self.opts.tmp_dir,suffix='rgFixMateOut') 
            self.delme.append(fixedbam)
            cl = ['java -Xmx%s' % self.opts.maxjheap,'-jar',self.opts.jar,'I=',opts.input,'O=',fixedbam] 
            # stitch the adjusted headers
            self.run(cl)
            sam = self.bamToSam(fixedbam)
            cl = ['cp',fixedbam,self.opts.output]
            self.run(cl)
        else:
            cl = ['java -Xmx%s' % self.opts.maxjheap,'-jar',self.opts.jar,'I=',opts.input,'O=',self.opts.output] 
            self.run(cl)
        
          
 
    def cleanUp(self):
        if self.tidy:
            for fname in self.delme:
                try:
                    os.unlink(fname)
                except:
                    pass
        self.tlog.close()
        s = open(self.tlogname,'r').readlines()
        print >> sys.stdout, s # for info
        os.unlink(self.tlogname)

    
if __name__ == '__main__':
    '''
    <command interpreter="python">
   rgPicardFixMate.py -i "$input_file" -n "$out_prefix" -o "$out_file" --tmp_dir "${__new_file_path__}"  
   -n "$newformat" -j "${GALAXY_DATA_INDEX_DIR}/shared/jars/FixMateInformation.jar"
   </command>
    '''
    op = optparse.OptionParser()
    op.add_option('-i', '--input', default=None)
    op.add_option('-o', '--output', default=None)
    op.add_option('-n', '--title', default="SortSamBam")
    op.add_option('--newformat', default='bam')
    op.add_option('-x', '--maxjheap', default='2g')
    op.add_option('--tmp_dir', default='/tmp')
    op.add_option('-j','--jar',default='')
    opts, args = op.parse_args()
    assert opts.input <> None
    assert os.path.isfile(opts.input) 
    assert os.path.isfile(opts.jar)
    try:
        os.makedirs(opts.tmp_dir)
    except:
        pass
    fm = fixMate(opts=opts)
    fm.doFix()
    fm.cleanUp()
    
    
