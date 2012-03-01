import os,string,sys,optparse,shutil,tempfile 
from subprocess import Popen 
from rgutils import galhtmlprefix,galhtmlpostfix,galhtmlattr,timenow, getFileString 
progname = os.path.split(sys.argv[0])[1]

""" 
insert size metrics 
Copyright ross lazarus oct 31 (happy halloween!) 
All rights reserved 
Released for rgenetics under the LGPL 
"""



class insertSize():
    """
    classy!
    """
    
    def __init__(self,opts=None,cl=[],tidy=False):
        """
        """
        self.ourname = 'rgPicardInsertSize'
        self.opts = opts
        self.tidy = tidy
        self.cl = ' '.join(cl) # ready for the htmlfile output
        self.delme = []
        killme = string.punctuation + string.whitespace
        trantab = string.maketrans(killme,'_'*len(killme))
        self.title = self.opts.title.translate(trantab)
        self.tlogname = os.path.join(self.opts.outdir,'rgInsertSizeMetrics.txt')
        self.tlog = open(self.tlogname,'w')
        self.isPDF = 'InsertSizeHist.pdf'
        self.info = '%s on %s at %s' % (self.ourname,self.title,timenow())
    

    def run(self,cl=None,redir=True):
        assert cl <> None
        fd,templog = tempfile.mkstemp(dir=self.opts.outdir,suffix='rgInsertSizeRun.txt')
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
        
    def bamToSam(self,infile=None):
        """
        use samtools view to convert bam to sam
        """
        fd,tempsam = tempfile.mkstemp(dir=self.opts.outdir,suffix='rgSortBamTemp.sam')
        cl = ['samtools view -h -o',tempsam,infile]
        self.run(cl)
        return tempsam

    def samToBam(self,infile=None):
        """
        use samtools view to convert sam to bam
        """
        fd,tempbam = tempfile.mkstemp(dir=self.opts.outdir,suffix='rgSortBamTemp.sam')
        cl = ['samtools view -h -b -S -o',tempbam,infile]
        self.run(cl)
        return tempbam
    
    def makehtml(self):
        """
        write the report as html
        """
        logdat = open(self.tlogname,'r').readlines()
        res = []
        res.append(galhtmlprefix % progname)
        res.append(galhtmlattr % (progname,timenow()))
        res.append('<b>Your job produced the following outputs - check here for a record of what was done and any unexpected events</b><hr/>')
        imghref = '%s.jpg' % os.path.splitext(self.isPDF)[0] # removes .pdf
        res.append('<table cellpadding="10"><tr><td>\n')
        res.append('<a href="%s"><img src="%s" alt="%s" hspace="10" align="middle"></a>\n' % (self.isPDF,imghref,imghref))
        res.append('</tr><td></table>\n')
        try:
            flist = os.listdir(self.opts.outdir)
        except:
            flist = []
        if len(flist) > 0: # we should clean everything up - picard doesn't tell us what it did in cleansam unfortunately
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
        res.append('orchestrated by the Galaxy rgInsertSize wrapper and this command line from the Galaxy form:<br/>\n%s' % (self.cl))
        res.append(galhtmlpostfix)
        f = open(self.opts.htmlout,'w')
        f.write('\n'.join(res))
        f.close()
  
        
    def doIS(self):
        """
        """
        hout = 'InsertSizeOut.txt'
        histpdf = self.isPDF
        cl = ['java -Xmx%s' % self.opts.maxjheap,'-jar',self.opts.jar,'I=',opts.input,'O=',hout,'HISTOGRAM_FILE=',histpdf]
        if self.opts.taillimit <> '0':
            cl.append('TAIL_LIMIT=%s' % self.opts.taillimit)
        if self.opts.histwidth <> '0':
            cl.append('HISTOGRAM_WIDTH=%s' % self.opts.histwidth)
        if float(self.opts.minpct) > 0.0:
            cl.append('MINIMUM_PCT=%s' %self.opts.minpct)
        self.run(cl)
        cl = ['mogrify', '-format jpg -resize x400 %s' % (os.path.join(self.opts.outdir,self.isPDF))]
        self.run(cl)        
        self.tlog.close()
           
 
    def cleanup(self):
        if self.tidy:
            for fname in self.delme:
                try:
                    os.unlink(fname)
                except:
                    pass
        print >> sys.stdout, self.info # for info

    
if __name__ == '__main__':
    '''
    <command interpreter="python">
   rgPicardInsertSize.py -i "$input_file" -n "$out_prefix" -o "$out_file" --tmp_dir "${__new_file_path__}"
   -l "$tailLimit" -w "$histWidth" -p "$minPct" -n "$newformat"
   -j "${GALAXY_DATA_INDEX_DIR}/shared/jars/CollectInsertSizeMetrics.jar" -d "$html_file.files_path" -t "$html_file"
  </command>
    '''
    op = optparse.OptionParser()
    op.add_option('-i', '--input', default=None)
    op.add_option('-n', '--title', default="Insert size metrics")
    op.add_option('-l', '--taillimit', default="0")
    op.add_option('-w', '--histwidth', default="0")
    op.add_option('-p', '--minpct', default="0.01")
    op.add_option('-t', '--htmlout', default="")
    op.add_option('-d', '--outdir', default="")
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
    try:
        os.makedirs(opts.outdir)
    except:
        pass
    x = insertSize(opts=opts,cl=sys.argv)
    x.doIS()
    x.makehtml()
    
