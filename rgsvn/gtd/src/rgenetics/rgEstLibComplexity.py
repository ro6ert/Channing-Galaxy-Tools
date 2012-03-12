import os,string,sys,optparse,shutil,tempfile 
from subprocess import Popen 
from rgutils import galhtmlprefix,galhtmlpostfix,galhtmlattr,timenow, getFileString 
progname = os.path.split(sys.argv[0])[1]

""" 
Estimate library complexity
Copyright ross lazarus oct 31 (happy halloween!) 
All rights reserved released for rgenetics under the LGPL 
"""



class libraryComplexity():
    """
    classy!
    """
    
    def __init__(self,opts=None,cl=[]):
        """
        """
        self.opts = opts
        self.cl = ' '.join(cl) # ready for the htmlfile output
        self.delme = []
        killme = string.punctuation + string.whitespace
        trantab = string.maketrans(killme,'_'*len(killme))
        self.title = self.opts.title.translate(trantab)
        self.tlogname = os.path.join(self.opts.outdir,'rgestimatelibcomp.txt')
        self.tlog = open(self.tlogname,'w')
        self.outtxt = 'estlibcompout.txt'
    

    def run(self,cl=None,redir=True):
        assert cl <> None
        fd,templog = tempfile.mkstemp(dir=self.opts.outdir,suffix='rgEstLibComplexity.txt')
        tlf = open(templog,'w')
        if redir:
           process = Popen(' '.join(cl), shell=True, stderr=tlf, stdout=tlf, cwd=self.opts.outdir)
        else:
           process = Popen(' '.join(cl), shell=True, cwd=self.opts.outdir)
        rval = process.wait()
        tlf.close()
        tlogs = ''.join(open(templog,'r').readlines())
        if len(tlogs) > 1:
            s = '## executing %s returned status %d and log (stdout/stderr) records: \n%s\n' % (' '.join(cl),rval,tlogs)
        else:
            s = '## executing %s returned status %d. Nothing appeared on stderr/stdout\n' % (' '.join(cl),rval)
        os.unlink(templog) # always
        self.tlog.write(s)        
    
    def makehtml(self):
        """
        write the report as html
        """
        logdat = open(self.tlogname,'r').readlines()
        res = []
        res.append(galhtmlprefix % progname)
        res.append(galhtmlattr % (progname,timenow()))
        res.append('<b>Your job produced the following outputs - check here for a record of what was done and any unexpected events</b><hr/>')
        try:
            flist = os.listdir(self.opts.outdir)
        except:
            flist = []
        if len(flist) > 0: # show what's left
            flist = [x for x in flist if not (x.startswith('.') or x == 'None')]
            tlist = [(os.path.getmtime(os.path.join(self.opts.outdir,x)),x) for x in flist]
            tlist.sort()
            flist = [x[1] for x in tlist]
            res.append('<div><b>Output files.</b><hr/>\n')
            res.append('<table>\n')
            for i,f in enumerate(flist):
                fn = os.path.split(f)[-1]
                fs = getFileString(fn,self.opts.outdir)
                res.append('<tr><td><a href="%s">%s</a></td></tr>\n' % (fn,fs))
            res.append('</table></div>\n')
        res.append('<b>Log of activity</b><hr/>\n')
        res.append('\n%s' % '<br/>'.join(logdat))
        res.append('<hr/>Note: The freely available <a href="http://picard.sourceforge.net/command-line-overview.shtml">Picard software</a> \n')
        res.append('generated all outputs reported here. These third party tools were')
        res.append('orchestrated by the Galaxy rgEstLibComplexity wrapper and this command line from the Galaxy form:<br/>\n%s' % (self.cl))
        res.append(galhtmlpostfix)
        f = open(self.opts.htmlout,'w')
        f.write('\n'.join(res))
        f.close()
  
        
    def doLC(self):
        """
        """

        cl = ['java -Xmx%s' % self.opts.maxjheap,'-jar',self.opts.jar,'I=',opts.input,'O=',self.outtxt]
        if float(self.opts.minid) > 0:
            cl.append('MIN_IDENTICAL_BASES=%s' % self.opts.minid)
        if float(self.opts.maxdiff) > 0.0:
            cl.append('MAX_DIFF_RATE=%s' % self.opts.maxdiff)
        if float(self.opts.minmeanq) > 0:
            cl.append('MIN_MEAN_QUALITY=%s' %self.opts.minmeanq)
        if self.opts.readregex > '':
            cl.append('READ_NAME_REGEX="%s"' %self.opts.readregex)
        if float(self.opts.optdupedist) > 0:
            cl.append('OPTICAL_DUPLICATE_PIXEL_DISTANCE=%s' %self.opts.optdupedist)
        self.run(cl)
        self.tlog.close()
           
if __name__ == '__main__':
    '''
    <command interpreter="python">
   rgEstLibComplexity.py -i "$input_file" -n "$out_prefix" --tmp_dir "${__new_file_path__}" --minid "$minID"
   --maxdiff "$maxDiff" --minmeanq "$minMeanQ" --readregex "$readRegex" --optdupedist "$optDupeDist"
   -j "${GALAXY_DATA_INDEX_DIR}/shared/jars/EstimateLibraryComplexity.jar" -d "$html_file.files_path" -t "$html_file"
   </command>

    '''
    op = optparse.OptionParser()
    op.add_option('-i', '--input', default=None)
    op.add_option('-n', '--title', default="Estimate Library Complexity")
    op.add_option('--minid', default="5")
    op.add_option('--maxdiff', default="0.03")
    op.add_option('--minmeanq', default="20")
    op.add_option('--readregex', default="[a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*")
    op.add_option('--optdupedist', default="100")
    op.add_option('-t', '--htmlout', default="")
    op.add_option('-d', '--outdir', default="")
    op.add_option('-x', '--maxjheap', default='2g')
    op.add_option('--tmp_dir', default='/tmp')
    op.add_option('-j','--jar',default='')
    opts, args = op.parse_args()
    assert opts.input <> None
    assert os.path.isfile(opts.input)
    assert os.path.isfile(opts.jar)
    if not os.path.exists(opts.tmp_dir): os.makedirs(opts.tmp_dir)
    if not os.path.exists(opts.outdir): os.makedirs(opts.outdir)

    tool = libraryComplexity(opts=opts,cl=sys.argv)
    tool.doLC()
    tool.makehtml()
    
