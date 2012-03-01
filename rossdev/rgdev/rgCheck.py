#!/usr/bin/env python
# test for current rgenetics tools
from rgutils import galhtmlprefix,plinke,RRun,runPlink 
# RRun(rcmd=[],outdir=None,title='myR')
# runPlink(bfn='bar',ofn='foo',logf='x.log',plinktasks=[],cd='./',vclbase = [])

def checkR(myTitle='Rgenetics installation test',outDir='./'):
    """ report missing packages or X11 for R """
    rexpLibs=['ArrayExpress','lumi','limma','AffyExpress',"affyQCReport",
    "GenABEL","simpleaffy","snpMatrix", "affyPLM",'haplo.stats','GEOquery',
    'arrayQualityMetrics']
    rgenLibs=['hexbin',]
    rchecks = ['X11()']
    res = []
    missinglibs = []
    reqLibs = rexpLibs + rgenLibs
    for thisLib in reqLibs: # ensure all R packages are available
        rlog,flist=RRun(['library(%s)' % thisLib,],outdir=outDir,title=myTitle)
        if rlog[0][:19] == 'Nonzero exit code =': # fail
            missinglibs.append(thisLib)
            res += rlog
    if len(missinglibs) > 0:
        res += ['The R package %s cannot be loaded - please install' % x for x in missinglibs] 
    rlog,flist = RRun(rchecks,outdir=outDir,title=myTitle)
    if rlog[0][:19] == 'Nonzero exit code =': # fail
            res.append('X11() failed in R - all graphical (png/jpg/pdf) outputs will fail')
            res.append('Was X11 configured at R compilation; is the virtual frame buffer Xvfb available for headless nodes?')
    return res
                
def checkUtils(myTitle="Rgenetics installation test",outDir='./'):
    """ check mogrify etc """
    utes = ['mogrify','pdfnup','convert','smartpca','java','pdfjoin','plink']
    res = []
    fplog,plog = tempfile.mkstemp()
    res.append('## Rgenetics: http://rgenetics.org Galaxy Tools Check for required utilities')
    for task in utes: # each is a list
        sto = file(plog,'w')
        x = subprocess.Popen(task,shell=True,cwd=cd,stdout=sto,stderr=sto)
        retval = x.wait()
        sto.close()
        try:
            lplog = file(plog,'r').read()
            os.unlink(plog) # no longer needed
        except:
            res.append('### %s Strange - no std out from %s' % (timenow(),task)
        if retval <> 0:
            res.append('%s not available - please install on nodes' % task)
        res += lplog
    return res
