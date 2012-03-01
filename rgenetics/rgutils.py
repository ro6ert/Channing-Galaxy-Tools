# utilities for rgenetics
#
# copyright 2009 ross lazarus
# released under the LGPL
#

import subprocess, os, sys, time, tempfile, string, plinkbinJZ, subprocess
import datetime

galhtmlprefix = """<?xml version="1.0" encoding="utf-8" ?>
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
<head>
<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
<meta name="generator" content="Galaxy %s tool output - see http://g2.trac.bx.psu.edu/" />
<title></title>
<link rel="stylesheet" href="/static/style/base.css" type="text/css" />
</head>
<body>
<div class="document">
"""
galhtmlattr = """<b><a href="http://rgenetics.org">Galaxy Rgenetics</a> tool output %s run at %s</b><br/>"""
galhtmlpostfix = """</div></body></html>\n"""

plinke = 'plink' # changed jan 2010 - all exes must be on path
rexe = 'R'       # to avoid cluster/platform dependencies
smartpca = 'smartpca.perl'

def timenow():
    """return current time as a string
    """
    return time.strftime('%d/%m/%Y %H:%M:%S', time.localtime(time.time()))

def timestamp():
    return datetime.datetime.now().strftime('%Y%m%d%H%M%S')

def fail( message ):
    print >> sys.stderr, message
    return -1

def whereis(program):
    for path in os.environ.get('PATH', '').split(':'):
        if os.path.exists(os.path.join(path, program)) and \
           not os.path.isdir(os.path.join(path, program)):
            return os.path.join(path, program)
    return None

def getFileString(fpath, outpath):
    """
    format a nice file size string
    """
    size = ''
    fp = os.path.join(outpath, fpath)
    s = '? ?'
    if os.path.isfile(fp):
        n = float(os.path.getsize(fp))
        if n > 2**20:
            size = ' (%1.1f MB)' % (n/2**20)
        elif n > 2**10:
            size = ' (%1.1f KB)' % (n/2**10)
        elif n > 0:
            size = ' (%d B)' % (int(n))
        s = '%s %s' % (fpath, size) 
    return s


class PicardBase():
    """
    simple base class with some utilities for Picard
    """
    
    def __init__(self, opts=None):
        """ common stuff needed at init for a picard tool
        """
        assert opts <> None, 'PicardBase needs opts at init'
        self.opts = opts
        assert self.opts.jar <> None, 'PicardBase needs opts to have a jar element'
        assert os.path.isfile(opts.input), 'PicardBase needs opts.input (%s supplied) that can be opened' % opts.input
        self.jar = self.opts.jar
        self.picname = os.path.splitext(os.path.basename(self.jar))[0] # get rid of .jar
        self.version = '0.001'
        killme = string.punctuation + string.whitespace
        self.trantab = string.maketrans(killme,'_'*len(killme))
        self.title = self.opts.title.translate(self.trantab)
        self.inputfile = opts.input
        try:
            os.makedirs(self.opts.outdir)
        except:
            pass
        try:
            os.makedirs(self.opts.tempdir)
        except:
            pass
        ifd, tempout = tempfile.mkstemp(prefix=self.picname,dir=self.opts.tmpdir)
        log_fd, self.log_filename = tempfile.mkstemp(dir=opts.tmpdir,suffix='log', prefix=self.picname)
        self.clparams = {} # key value pairs for constructing variable parts of tool command line
        self.delme = [] # list of files to destroy
       
    def runCL(self,cl=None,output_dir=None):
        """ construct and run a command line
        """
        assert cl <> None, 'PicardBase runCL needs a command line as cl'
        fd,templog = tempfile.mkstemp(dir=output_dir,suffix='rgtempRun.txt')
        tlf = open(templog,'w')
        if output_dir <> None:
           process = subprocess.Popen(' '.join(cl), shell=True, stderr=tlf, stdout=tlf, cwd=output_dir)
        else:
           process = subprocess.Popen(' '.join(cl), shell=True)
        rval = process.wait()
        tlf.close()
        tlogs = ''.join(open(templog,'r').readlines())
        if len(tlogs) > 1:
            s = '## executing %s returned status %d and log (stdout/stderr) records: \n%s\n' % (' '.join(cl),rval,tlogs)
        else:
            s = '## executing %s returned status %d. Nothing appeared on stderr/stdout\n' % (' '.join(cl),rval)
        os.unlink(templog) # always
        return s
    
    def runPic(self):
        """
        self.clparams should be a dict of key value pairs to write to the command line
        params = {'I=':self.opts.input}
        params['O='] = tempout
        """
        cll = ['java -Xmx%s' % self.opts.maxjheap,' -jar %s' % self.jar]
        cll += [''.join(it) for it in list(self.clparams.items())]
        self.cl = ' '.join(cll)
        s = self.runCL(cl=cll,output_dir=self.opts.outdir)
        lf = open(self.log_filename,'a')
        lf.write(s)
        lf.write('\n')
        lf.close()

    def samToBam(self,infile=None,outdir=None):
        """
        use samtools view to convert sam to bam
        """
        fd,tempbam = tempfile.mkstemp(dir=outdir,suffix='rgSortBamTemp.bam')
        cl = ['samtools view -h -b -S -o ',tempbam,infile]
        tlog = self.runCL(cl,outdir)
        return tlog,tempbam
    
    def bamToSam(self,infile=None,outdir=None):
        """
        use samtools view to convert bam to sam
        """
        fd,tempsam = tempfile.mkstemp(dir=outdir,suffix='rgutilsTemp.bam')
        cl = ['samtools view -h -o ',tempsam,infile]
        tlog = self.runCL(cl,outdir)
        return tlog,tempsam

    def sortSam(self, infile=None,outfile=None,outdir=None):
        """
        """
        print '## sortSam got infile=%s,outfile=%s,outdir=%s' % (infile,outfile,outdir)
        cl = ['samtools sort',infile,outfile]
        tlog = self.runCL(cl,outdir)
        return tlog

    def cleanup(self):
        for fname in self.delme:
            try:
                os.unlink(fname)
            except:
                pass
                    
    def prettyPicout(self,picout=None,transpose=True,maxrows=100):
        """organize picard outpouts into a report html page
        """
        res = []
        try:
            r = open(picout,'r').readlines()
        except:
            r = []        
        if len(r) > 0:
            res.append('<b>Picard on line resources</b><ul>\n')
            res.append('<li><a href="http://picard.sourceforge.net/index.shtml">Click here for Picard Documentation</a></li>\n')
            res.append('<li><a href="http://picard.sourceforge.net/picard-metric-definitions.shtml">Click here for Picard Metrics definitions</a></li></ul><hr/>\n')
            if transpose:
                res.append('<b>Picard output (transposed to make it easier to see)</b><hr/>\n')       
            else:
                res.append('<b>Picard output</b><hr/>\n')  
            res.append('<table cellpadding="3" >\n')
            dat = []
            heads = []
            lastr = len(r) - 1
            # special case for estimate library complexity hist
            thist = False
            for i,row in enumerate(r):
                if row.strip() > '':
                    srow = row.split('\t')
                    if row.startswith('#'):
                        heads.append(row.strip()) # want strings
                    else:
                        dat.append(srow) # want lists
                    if row.startswith('## HISTOGRAM'):
                        thist = True
                if row.strip() == '' or i == lastr: # last line or blank means write a segment
                    if len(heads) > 0:
                        hres = ['<tr class="d%d"><td>%s</td></tr>' % (i % 2,x) for i,x in enumerate(heads)]
                        res += hres
                        heads = []
                    if len(dat) > 0:  
                        if transpose and not thist:
                            tdat = map(None,*dat) # transpose an arbitrary list of lists
                            tdat = ['<tr class="d%d"><td>%s</td><td>%s</td></tr>\n' % ((i+len(heads)) % 2,x[0],x[1]) for i,x in enumerate(tdat) if i < maxrows] 
                            missing = len(tdat) - maxrows
                            if missing > 0:
                               tdat.append('<tr><td colspan="2">...WARNING: %d rows deleted..see raw file %s for entire output</td></tr>' % (missing,os.path.basename(picout)))
                        else:
                            if thist:
                                dat = ['\t'.join(x).strip() for x in dat] # back to strings :(
                            tdat = ['<tr class="d%d"><td>%s</td></tr>\n' % ((i+len(heads)) % 2,x) for i,x in enumerate(dat) if i < maxrows] 
                            missing = len(tdat) - maxrows
                            if missing > 0:      
                                tdat.append('<tr><td>...WARNING: %d rows deleted..see raw file %s for entire output</td></tr>' % (missing,os.path.basename(picout)))
                        res += tdat
                        dat = []
            res.append('</table>\n')   
        else:
            res.append('<b>No Picard output found - please consult the Picard log below for an explanation</b><br/>\n')
        return res

    def fixPicardOutputs(self,tempout=None,output_dir=None,log_file=None,html_output=None,progname=None,transpose=True,maxloglines=100):
        """
        picard produces long hard to read tab header files
        make them available but present them transposed for readability
        """
        self.cleanup() # remove temp files stored in delme
        self.maxloglines = maxloglines
        rstyle="""<style type="text/css">
        tr.d0 td {background-color: oldlace; color: black;}
        tr.d1 td {background-color: aliceblue; color: black;}
        </style>"""    
        res = [rstyle,]
        res.append(galhtmlprefix % progname)   
        res.append(galhtmlattr % (progname,timenow()))
        flist = [x for x in os.listdir(output_dir) if not x.startswith('.')] 
        pdflist = [x for x in flist if os.path.splitext(x)[-1].lower() == '.pdf']
        if len(pdflist) > 0: # assumes all pdfs come with thumbnail .jpgs
            for p in pdflist:
                imghref = '%s.jpg' % os.path.splitext(p)[0] # removes .pdf
                res.append('<table cellpadding="10"><tr><td>\n')
                res.append('<a href="%s"><img src="%s" title="Click image preview for a print quality PDF version" hspace="10" align="middle"></a>\n' % (p,imghref)) 
                res.append('</tr></td></table>\n')   
        if len(flist) > 0:
            res.append('<b>Running this Galaxy tool produced the following output files (click the filename to view/download a copy).</b><hr/>')
            res.append('<table>\n')
            for i,f in enumerate(flist):
                fn = os.path.split(f)[-1]
                res.append('<tr><td><a href="%s">%s</a></td></tr>\n' % (fn,fn))
            res.append('</table><p/>\n') 
        pres = self.prettyPicout(picout=tempout,transpose=transpose,maxrows=self.maxloglines)
        res += pres
        l = open(self.log_filename,'r').readlines()
        llen = len(l)
        if llen > 0: 
            res.append('<b>Picard log</b><hr/>\n') 
            rlog = ['<pre>',]
            if llen > self.maxloglines:
                rlog += l[:self.maxloglines]                
                rlog.append('\n<b>## WARNING - %d log lines truncated - %s contains entire output' % (llen - self.maxloglines,self.log_filename))
            else:
                rlog += l
            rlog.append('</pre>')
            res += rlog
        else:
            res.append("### Odd, Picard left no log file %s - must have really barfed badly?\n" % self.log_filename)
        res.append('<hr/>The freely available <a href="http://picard.sourceforge.net/command-line-overview.shtml">Picard software</a> \n') 
        res.append( 'generated all outputs reported here, using this command line:<br/>\n<pre>%s</pre>\n' % self.cl)   
        res.append(galhtmlpostfix) 
        outf = open(html_output,'w')
        outf.write(''.join(res))   
        outf.write('\n')
        outf.close()

    def bedToPicInterval(self, infile=None, outf=None):
        """
        deprecated - kept because it's a way of recovering informative sam style headers for bed files
        we now use a simple cut and paste of the incoming bam header to head of the bait/target files
        because otherwise we lose our minds:

        Aaaarrrrrggggghhhh - the bait/target header that picard wants 
        (incidentally, how does this help anyone? Why not just read a bed format file?)
        MUST have the same number of sequences as the incoming bam/sam
        Bam header (samtools view -H) pasted to bait/target works fine to fix the message below:

[Wed Apr 06 11:53:46 EDT 2011] net.sf.picard.analysis.directed.CalculateHsMetrics 
BAIT_INTERVALS=/share/shared/galaxy/database/job_working_directory/2165/dataset_2895_files/rgPicardHsMetrics.bait 
TARGET_INTERVALS=/share/shared/galaxy/database/job_working_directory/2165/dataset_2895_files/rgPicardHsMetrics.bait 
INPUT=/share/shared/galaxy/database/files/002/dataset_2862.dat 
OUTPUT=/share/shared/galaxy/database/job_working_directory/2165/dataset_2895_files/rgPicardHsMetrics.out TMP_DIR=/share/shared/galaxy/database/tmp 
VALIDATION_STRINGENCY=LENIENT    VERBOSITY=INFO QUIET=false COMPRESSION_LEVEL=5 MAX_RECORDS_IN_RAM=500000 CREATE_INDEX=false CREATE_MD5_FILE=false
[Wed Apr 06 11:53:52 EDT 2011] net.sf.picard.analysis.directed.CalculateHsMetrics done.
Runtime.totalMemory()=471531520
Exception in thread "main" net.sf.samtools.util.SequenceUtil$SequenceListsDifferException: Sequence dictionaries are not the same size (25, 24)
	at net.sf.samtools.util.SequenceUtil.assertSequenceListsEqual(SequenceUtil.java:93)
	at net.sf.samtools.util.SequenceUtil.assertSequenceDictionariesEqual(SequenceUtil.java:135)
	at net.sf.picard.analysis.directed.CalculateHsMetrics.doWork(CalculateHsMetrics.java:87)
	at net.sf.picard.cmdline.CommandLineProgram.instanceMain(CommandLineProgram.java:157)
	at net.sf.picard.analysis.directed.CalculateHsMetrics.main(CalculateHsMetrics.java:68)

        Picard tools requiring targets want
        a sam style header which incidentally, MUST be sorted in natural order - not lexicographic order:

        @SQ     SN:chrM LN:16571
        @SQ     SN:chr1 LN:247249719
        @SQ     SN:chr2 LN:242951149
        @SQ     SN:chr3 LN:199501827
        @SQ     SN:chr4 LN:191273063
        added to the start of what looks like a bed style file
        chr1    67052400        67052451        -       CCDS635.1_cds_0_0_chr1_67052401_r
        chr1    67060631        67060788        -       CCDS635.1_cds_1_0_chr1_67060632_r
        chr1    67065090        67065317        -       CCDS635.1_cds_2_0_chr1_67065091_r
        chr1    67066082        67066181        -       CCDS635.1_cds_3_0_chr1_67066083_r


        see http://genome.ucsc.edu/FAQ/FAQtracks.html#tracks1
        we need to add 1 to start coordinates on the way through - but length calculations are easier
        # bedToPicard.py
        # ross lazarus October 2010
        # LGPL
        # for Rgenetics
        """

        def getFlen(bedfname=None):
            """
            find all features in a BED file and sum their lengths
            """
            features = {}
            infile = open(bedfname,'r')
            for i,row in enumerate(infile):
                if row[0] == '@': # shouldn't happen given a bed file!
                    print >> sys.stderr, 'row %d=%s - should NOT start with @!' % (i,row)
                    sys.exit(1)
                row = row.strip()
                if len(row) > 0:
                    srow = row.split('\t')
                    f = srow[0]
                    spos = srow[1] # zero based from UCSC so no need to add 1 - eg 0-100 is 100 bases numbered 0-99 (!)
                    epos = srow[2] # see http://genome.ucsc.edu/FAQ/FAQtracks.html#tracks1
                    flen = int(epos) - int(spos)
                    features.setdefault(f,0)
                    features[f] += flen
            infile.close()
            return features

        def keynat(string):
            '''
            borrowed from http://code.activestate.com/recipes/285264-natural-string-sorting/
            A natural sort helper function for sort() and sorted()
            without using regular expressions or exceptions.

            >>> items = ('Z', 'a', '10th', '1st', '9')
            >>> sorted(items)
            ['10th', '1st', '9', 'Z', 'a']
            >>> sorted(items, key=keynat)
            ['1st', '9', '10th', 'a', 'Z']
            '''
            it = type(1)
            r = []
            for c in string:
                if c.isdigit():
                    d = int(c)
                    if r and type( r[-1] ) == it:
                        r[-1] = r[-1] * 10 + d
                    else:
                        r.append(d)
                else:
                    r.append(c.lower())
            return r

        def writePic(outfname=None,bedfname=None):
            """
            collect header info and rewrite bed with header for picard
            """
            featlen = getFlen(bedfname=bedfname)
            outf = open(outfname,'w')
            infile = open(bedfname,'r') # already tested in getFlen
            k = featlen.keys()
            fk = sorted(k, key=keynat)
            header = ['@SQ\tSN:%s\tLN:%d' % (x,featlen[x]) for x in fk]
            outf.write('\n'.join(header))
            outf.write('\n')
            for row in infile:
                row = row.strip()
                if len(row) > 0: # convert zero based start coordinate to 1 based
                    srow = row.split('\t')
                    srow[1] = '%d' % (int(srow[1])+1) # see http://genome.ucsc.edu/FAQ/FAQtracks.html#tracks1
                    outf.write('\t'.join(srow))
                    outf.write('\n')
            outf.close()
            infile.close()
            return featlen



        # bedToPicInterval starts here
        assert os.path.isfile(infile),'## bedtoPicInterval error - cannot open supplied input bed file %s' % infile
        if outf == None:
            fd,outf = tempfile.mkstemp(prefix='rgPicardHsMetrics')
        writePic(outfname=outf,bedfname=infile)
        return outf

def RRun(rcmd=[],outdir=None,title='myR',tidy=True):
    """
    run an r script, lines in rcmd,
    in a temporary directory
    move everything, r script and all back to outdir which will be an html file


      # test
      RRun(rcmd=['print("hello cruel world")','q()'],title='test')
    echo "a <- c(5, 5); b <- c(0.5, 0.5)" | cat - RScript.R | R --slave \ --vanilla
    suggested by http://tolstoy.newcastle.edu.au/R/devel/05/09/2448.html
    """
    killme = string.punctuation + string.whitespace
    trantab = string.maketrans(killme,'_'*len(killme))
    title = title.translate(trantab)
    rlog = []
    tempout=False
    rname = '%s.R' % title
    stoname = '%s.R.log' % title
    cwd = os.getcwd()
    if outdir: # want a specific path
        try:
            os.makedirs(outdir) # might not be there yet...
        except:
            pass
        os.chdir(outdir)
    if type(rcmd) == type([]):
        script = '\n'.join(rcmd)
    else: # string
        script = rcmd
    sto = file(stoname,'w')
    rscript = file(rname,'w')
    rscript.write(script)
    rscript.write('\n#R script autogenerated by rgenetics/rgutils.py on %s\n' % timenow())
    rscript.close()
    vcl = '%s --slave --vanilla < %s' %  (rexe,rname)
    if outdir:
        x = subprocess.Popen(vcl,shell=True,stderr=sto,stdout=sto,cwd=outdir)
    else:
        x = subprocess.Popen(vcl,shell=True,stderr=sto,stdout=sto)
    retval = x.wait()
    sto.close()
    rlog = file(stoname,'r').readlines()
    if retval <> 0:
        rlog.insert(0,'Nonzero exit code = %d' % retval) # indicate failure
    if outdir:
        flist = os.listdir(outdir)
    else:
        flist = os.listdir(os.getcwd())
    flist.sort
    flist = [(x,x) for x in flist]
    for i,x in enumerate(flist):
        if x == rname:
            flist[i] = (x,'R script for %s' % title)
        elif x == stoname:
            flist[i] = (x,'R log for %s' % title)
    if outdir:
        os.chdir(cwd)
    return rlog,flist # for html layout

def runPlink(bfn='bar',ofn='foo',logf=None,plinktasks=[],cd='./',vclbase = []):
    """run a series of plink tasks and append log results to stdout
    vcl has a list of parameters for the spawnv
    common settings can all go in the vclbase list and are added to each plinktask
    """
    # root for all
    fplog,plog = tempfile.mkstemp()
    if type(logf) == type('  '): # open otherwise assume is file - ugh I'm in a hurry
    	mylog = file(logf,'a+')
    else:
        mylog = logf
    mylog.write('## Rgenetics: http://rgenetics.org Galaxy Tools rgQC.py Plink runner\n')
    for task in plinktasks: # each is a list
        vcl = vclbase + task
        sto = file(plog,'w')
        x = subprocess.Popen(' '.join(vcl),shell=True,stdout=sto,stderr=sto,cwd=cd)
        retval = x.wait()
        sto.close()
        try:
            lplog = file(plog,'r').read()
            mylog.write(lplog)
            os.unlink(plog) # no longer needed
        except:
            mylog.write('### %s Strange - no std out from plink when running command line\n%s' % (timenow(),' '.join(vcl)))

def pruneLD(plinktasks=[],cd='./',vclbase = []):
    """
    plink blathers when doing pruning - ignore
    Linkage disequilibrium based SNP pruning
    if a million snps in 3 billion base pairs, have mean 3k spacing
    assume 40-60k of ld in ceu, a window of 120k width is about 40 snps
    so lots more is perhaps less efficient - each window computational cost is
    ON^2 unless the code is smart enough to avoid unecessary computation where
    allele frequencies make it impossible to see ld > the r^2 cutoff threshold
    So, do a window and move forward 20?
    The fine Plink docs at http://pngu.mgh.harvard.edu/~purcell/plink/summary.shtml#prune
    reproduced below

Sometimes it is useful to generate a pruned subset of SNPs that are in approximate linkage equilibrium with each other. This can be achieved via two commands: 
--indep which prunes based on the variance inflation factor (VIF), which recursively removes SNPs within a sliding window; second, --indep-pairwise which is 
similar, except it is based only on pairwise genotypic correlation.

Hint The output of either of these commands is two lists of SNPs: those that are pruned out and those that are not. A separate command using the --extract or 
--exclude option is necessary to actually perform the pruning.

The VIF pruning routine is performed:
plink --file data --indep 50 5 2

will create files

     plink.prune.in
     plink.prune.out

Each is a simlpe list of SNP IDs; both these files can subsequently be specified as the argument for
a --extract or --exclude command.

The parameters for --indep are: window size in SNPs (e.g. 50), the number of SNPs to shift the
window at each step (e.g. 5), the VIF threshold. The VIF is 1/(1-R^2) where R^2 is the multiple correlation coefficient for a SNP being regressed on all other 
SNPs simultaneously. That is, this considers the correlations between SNPs but also between linear combinations of SNPs. A VIF of 10 is often taken to represent 
near collinearity problems in standard multiple regression analyses (i.e. implies R^2 of 0.9). A VIF of 1 would imply that the SNP is completely independent of 
all other SNPs. Practically, values between 1.5 and 2 should probably be used; particularly in small samples, if this threshold is too low and/or the window 
size is too large, too many SNPs may be removed.

The second procedure is performed:
plink --file data --indep-pairwise 50 5 0.5

This generates the same output files as the first version; the only difference is that a
simple pairwise threshold is used. The first two parameters (50 and 5) are the same as above (window size and step); the third parameter represents the r^2 
threshold. Note: this represents the pairwise SNP-SNP metric now, not the multiple correlation coefficient; also note, this is based on the genotypic 
correlation, i.e. it does not involve phasing.

To give a concrete example: the command above that specifies 50 5 0.5 would a) consider a
window of 50 SNPs, b) calculate LD between each pair of SNPs in the window, b) remove one of a pair of SNPs if the LD is greater than 0.5, c) shift the window 5 
SNPs forward and repeat the procedure.

To make a new, pruned file, then use something like (in this example, we also convert the
standard PED fileset to a binary one):
plink --file data --extract plink.prune.in --make-bed --out pruneddata
    """
    fplog,plog = tempfile.mkstemp()
    alog = []
    alog.append('## Rgenetics: http://rgenetics.org Galaxy Tools rgQC.py Plink pruneLD runner\n')
    for task in plinktasks: # each is a list
        vcl = vclbase + task
        sto = file(plog,'w')
        x = subprocess.Popen(' '.join(vcl),shell=True,stdout=sto,stderr=sto,cwd=cd)
        retval = x.wait()
        sto.close()
        try:
            lplog = file(plog,'r').readlines()
            lplog = [x for x in lplog if x.find('Pruning SNP') == -1]
            alog += lplog
            alog.append('\n')
            os.unlink(plog) # no longer needed
        except:
            alog.append('### %s Strange - no std out from plink when running command line\n%s\n' % (timenow(),' '.join(vcl)))
    return alog

def readMap(mapfile=None,allmarkers=False,rsdict={},c=None,spos=None,epos=None):
    """abstract out - keeps reappearing
    """
    mfile = open(mapfile, 'r')
    markers = []
    snpcols = {}
    snpIndex = 0 # in case empty or comment lines
    for rownum,row in enumerate(mfile):
        line = row.strip()
        if not line or line[0]=='#': continue
        chrom, snp, genpos, abspos = line.split()[:4] # just in case more cols
        try:
            abspos = int(abspos)
        except:
            abspos = 0 # stupid framingham data grumble grumble
        if allmarkers or rsdict.get(snp,None) or (chrom == c and (spos <= abspos <= epos)):
            markers.append((chrom,abspos,snp)) # decorate for sort into genomic
            snpcols[snp] = snpIndex # so we know which col to find genos for this marker
            snpIndex += 1
    markers.sort()
    rslist = [x[2] for x in markers] # drop decoration
    rsdict = dict(zip(rslist,rslist))
    mfile.close()
    return markers,snpcols,rslist,rsdict


