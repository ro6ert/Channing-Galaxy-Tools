import  sys, string, os, tempfile, shutil, subprocess, glob,time

# rpy removed sept 2009 rml
# replaced with direct RRun call to os.popen..
# eesh. ArrayExpress occasionally (eg E-MEXP-156) 
# returns variable number of different datatypes
# multiple chips = multipe affybatch/eset 
# ? need phenotype for each one
# should split multichip experiments? 
# ? make an rexset - an html derived type that can report presence or convert
# from the ur type - affybatch to eset or MAlist or whatever - or create these all
# when build a new one?
# converted to rpy in order to setup for variable output files

# first try at a BioC R script
# copyright 2008, ross lazarus 
# started August 6 2008
# licensed under the LGPL http://www.gnu.org/licenses/lgpl.html
# Not fit for human consumption
#
# Notes and code mostly cribbed from the BioC GEOquery vignette
#
# this is a non-interactive script intended as the
# executable for a Galaxy tool
# to create a new ExpressionSet and MAList from a GEO GDSxxxx id
# using the GEOquery package from BioC
#
# finally figured out how to finesse R's blathering output to stderr
# which makes the current Galaxy job runner think the job has failed
#
# http://www.bioconductor.org/workshops/2005/BioC2005/labs/lab01/index.html
# has good usecases.

raecode = """
doLocal<-function(eset,outprefix,tdir) 
{
    outf_pheno=paste(outprefix,'pheno',sep='.')
    outf_features=paste(outprefix,'features',sep='.')
    # need to get rid of :,;~| and other popular delimiters
    # these break , delimited parameter passing so must be dealt
    # whenever we save a file - no idea where you came from!
    rexphead = "[^(a-z,A-Z,0-9,'.')]" # allowed in names
    rexpval = "[^(a-z,A-Z,0-9,'.',':','/','_')]" # allowed in values and names
    phe = pData(eset)
    vm = varMetadata(eset)
    for ( i in c(1:length(names(phe))) )
       { phe[,i] = gsub(rexpval, "_", phe[,i]) }
    # remove all non letters/numbers
    newp = new("AnnotatedDataFrame",data=phe,varMetadata=vm)
    eset@phenoData = newp
    ppath = file.path(tdir,outf_pheno)
    write.table(phe,file=ppath,sep='\t',quote=F,row.names=F) 
    fpath = file.path(tdir,outf_features)
    fdd = eset@featureData@data
    write.table(fdd,file=fpath,sep='\t',quote=F,row.names=F)
    return(eset)
} # doLocal


mysaveAE = function(ae=NULL,AEId='None',tdir='.') {
 outf_eset=paste(AEId,'eset',sep='.')
 outf_affybatch=paste(AEId,'affybatch',sep='.')
 outf_nchannel=paste(AEId,'nchannel',sep='.')

 if (class(ae) == 'list') 
 # this is a right royal pain
 {
    for (i in 1:length(ae))
    {
        aepart = ae[[i]] 
        # retrieve list part
        klass = class(aepart)[1] 
        # eg "AffyBatch"
        print(paste('i=',i,'class(ae[[i]])=',klass))
        outprefix = paste(AEId,i,sep='_')
        aepart = doLocal(eset=aepart,outprefix=outprefix,tdir=tdir)   
        if (klass == 'AffyBatch') 
            {
            aef = file.path(tdir,paste(i,outf_affybatch,sep='_'))
            save(aepart,file=aef)        
            nae = rma(aepart) 
            # normalize
            naef=file.path(tdir,paste(i,outf_eset,sep='_'))
            save(nae,file=naef)  
            }
        if (klass == 'ExpressionSet' | klass == 'lumi')
            {
            aef = file.path(tdir,paste(i,outf_eset,sep='_'))
            save(aepart,file=aef) 

            }
        if (klass == 'NChannel')
            {
            ncf = file.path(tdir,paste(i,outf_nchannel,sep='_'))
            save(aepart,file=ncf)

            }
    } 
    # end for ae list of objects
 } # list
 else
 {
    klass = class(ae)[1]
    print(paste('class(ae)=',klass))
    outprefix = AEId
    ae = doLocal(eset=ae,outprefix=outprefix,tdir=tdir)   
    if (klass == 'AffyBatch') 
        {
        aef = file.path(tdir,outf_affybatch)
        save(ae,file=aef)        
        nae = rma(ae) # normalize
        naef = file.path(tdir,outf_eset)
        save(nae,file=naef)  
        }
    if (klass == 'ExpressionSet')
        {
        aef = file.path(tdir,outf_eset)
        save(ae,file=aef)        

        }
    if (klass == 'NChannel')
        {
        aef = file.path(tdir,outf_nchannel)
        save(ae,file=aef)

        }

 }

} # saveae - called for raw and processed separately March 2010
# need to know what was written to make galaxy datasets


rae = function(AEId='%s',tdir='%s')
{
# designed to write a series of whatever comes out of an arrayexpress id
# writes a list of incrementing filenames and a list of their
# types for each repeated element to the log file
# for the after-job hook to use to create new Galaxy datatypes
# this is part of the rexpression Galaxy toolkit
# copyright ross lazarus September 2008
# licensed under the lgpl
#

cat(paste('## tdir = ',tdir,'\n',sep=""))
library('ArrayExpress')
library('affy')
library('lumi')

cat(paste('### Array Express id =',AEId))

ae = getAE(AEId, tdir, type = 'full')
## Build a an ExpressionSet from the processed data
aecols = getcolproc(ae)
aeset = try(procset(ae,aecols[2]))
if (!inherits(aeset,"try-error")) {
     mysaveAE(ae=aeset,AEId=AEId,tdir=tdir)
     print(paste('## Saved processed data for',AEId,'\n'))
     } else {
       print(paste('### Warning - ArrayExpress unable to obtain processed data for',AEId))
     }
## Build a an ExpressionSet from the raw data
aeraw = try(magetab2bioc(files = ae))
if (!inherits(aeraw,"try-error")) {
     mysaveAE(ae=aeraw,AEId=AEId,tdir=tdir)
     print(paste('## Saved raw data for',AEId,'\n'))
     } else {
       print(paste('### Warning - ArrayExpress unable to obtain raw data for',AEId))
     }
} 
# end function rae


rae() # call with defaults
"""
# instantiate this with % (AEid,tdir)
# to run



galhtmlprefix = """<?xml version="1.0" encoding="utf-8" ?>
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
<head>
<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
<meta name="generator" content="Galaxy %s tool output - see http://galaxyproject.org/" />
<title></title>
<link rel="stylesheet" href="/static/style/base.css" type="text/css" />
</head>
<body>
<div class="document">
"""



galhtmlpostfix = """</div></body></html>"""
debug = False # useful for debugging RRun()

def timenow():
    """return current time as a string
    """
    return time.strftime('%d/%m/%Y %H:%M:%S', time.localtime(time.time()))



def whereis(program):
    for path in os.environ.get('PATH', '').split(':'):
        if os.path.exists(os.path.join(path, program)) and \
           not os.path.isdir(os.path.join(path, program)):
            return os.path.join(path, program)
    return None


def RRun(rcmd=[],outdir='.',title='myR',rexe=''):
    """
    run an r script, lines in rcmd,
    in a temporary directory
    move everything, r script and all back to outdir which will be an html file

    
      # test
      RRun(rcmd=['print("hello cruel world")','q()'],title='test')    
    """
    rlog = []
    assert os.path.isfile(rexe) 
    rname = '%s.R' % title
    stoname = '%s.R.log' % title 
    rfname = os.path.join(outdir,rname)
    try:
        os.makedirs(outdir) # might not be there yet...
    except:
        pass
    f = open(rfname,'w')
    if type(rcmd) == type([]):
        f.write('\n'.join(rcmd))
    else: # string
        f.write(rcmd)
    f.write('\n')
    f.close()
    stofname = os.path.join(outdir,stoname)
    sto = file(stofname,'w')
    vcl = [rexe,"--vanilla --slave", '<', rfname ]
    x = subprocess.Popen(' '.join(vcl),shell=True,stderr=sto,stdout=sto,cwd=outdir)
    retval = x.wait()
    sto.close()
    rlog = file(stofname,'r').readlines()
    flist = os.listdir(outdir)
    flist.sort
    flist = [(x,x) for x in flist] 
    for i,x in enumerate(flist):
        if x == rname:
            flist[i] = (x,'## R script for %s' % title)
        elif x == stoname:
            flist[i] = (x,'## R log for %s' % title)        
    return rlog,flist # for html layout
 

def main():
    """
<command interpreter="python">
reArrayExpress.py '$AEid'  '$title' '$htmlfile' '$htmlfile.extra_files_path' '${GALAXY_DATA_INDEX_DIR}/rg/bin/R'
</command>
    """
    nparm = 5
    appname = sys.argv[0]
    if (len(sys.argv) < nparm):
        print '%s needs %d parameters - given %d as %s' % (appname,nparm,len(sys.argv),';'.join(sys.argv))
        sys.exit(1)
        # short command line error
    appname = sys.argv[0]
    AEid = sys.argv[1].strip()
    title = sys.argv[2].strip()
    logfname = sys.argv[3].strip()
    outdir = sys.argv[4].strip()
    rexe = sys.argv[5].strip()
    try:
        os.makedirs(outdir) # may not exist yet
    except:
        pass
    replace = string.whitespace + string.punctuation
    ttab = string.maketrans(replace,'_'*len(replace))
    title = string.translate(title,ttab)
    tdir = tempfile.mkdtemp(prefix=AEid)
    logf = file(logfname,'w')
    logf.write(galhtmlprefix % appname)
    logf.write('<div># %s - part of the Rexpression Galaxy toolkit http://esphealth.org/trac/rgalaxy<br>\n' % (appname))
    logf.write('# Started at %s with parameters %s<br></div>\n' % (timenow(),' '.join(['"%s"' % x for x in sys.argv])))
    # AEId,outf_affybatch,outf_eSet,outf_nchannel,outf_pheno
    R = raecode % (AEid,outdir)
    rlog,flist = RRun(rcmd=R,outdir=outdir,title=title,rexe=rexe)
    flist.sort()
    logf.write('<div><table cellpadding="5" border="0">')
    for (fname,blurb) in flist:
        fpath = os.path.join(outdir,fname)
        size = ''
        if os.path.isfile(fpath):
            size = '(%1.1eMB)' % (float(os.path.getsize(fpath))/2**20)
        logf.write('<tr><td><a href="%s">%s %s</a></td></tr>' % (fname,blurb,size)) # or so
    logf.write('</table><hr/></div>')
    logf.write('\n<pre>\n %s' % ''.join(rlog)) 
    logf.write('\n</pre>')
    logf.write(galhtmlpostfix)
    logf.close()

if __name__ == "__main__":
    main()
