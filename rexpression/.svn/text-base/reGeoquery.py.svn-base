# reGeoquery.py
# started out with Rpy
# but clean solution with reusable R artefact
# that writes a generic R function and calls it
# at the end of the script with the form parameters
# where they should be in the call.
# copyright ross lazarus 2008
# released under the LGPL
# for the rexpression project
#
# sept 2009
# eesh. some metadata from geo and arrayexpress break the clumsy
# string handling - eg many have embedded delimeters like ','
# so our dumb parameter passing fails miserably
# what to do - improve parameter passing?
# uses BioC GEOquery
# deprecated: (note tricks to return output files encoded in the log file)
# now uses: reMakeNew_code.py for exec_after_process hook
# to create arbitrary rgenetics datasets as history outputs
#
# derived from arrayexpress code to return variable/unknown number
# of rexpression datasets with metadata
# these are created in reGeoquery_code.py
#
# updated jan 2009 to try to get genvar GSE6536
"""
# ross lazarus sept 2009
# problem getting a large human dataset split into parts with identical columns - so need to append the second one?
library(GEOquery)
getGEO('GSE8052')
# failed - sept 13 - ftp probem - wget worked fast - fixed on sept 14 go figure
# wget 'ftp://ftp.ncbi.nih.gov/pub/geo/DATA/SeriesMatrix/GSE8052/GSE8052_series_matrix-1.txt.gz'
# wget 'ftp://ftp.ncbi.nih.gov/pub/geo/DATA/SeriesMatrix/GSE8052/GSE8052_series_matrix-2.txt.gz'
#f1 = '/home/rerla/rcode/GSE8052_series_matrix-1.txt.gz'
#gse8052 = getGEO(file=f1)
#f2 = '/home/rerla/rcode/GSE8052_series_matrix-2.txt.gz'
#gse8052_2 = getGEO(file=f2)
gse8052_2 = gse8052[[2]]
gse8052 = gse8052[[1]]
a1 = gse8052@assayData
a2 = gse8052_2@assayData
e1 = a1$exprs
e2 = a2$exprs
anno = a1@annotation
gse8052@annotation = anno
newe = cbind(e1,ea2)
gse8052@assayData = assayDataNew("lockedEnvironment", exprs=newe)
p1 = pData(gse8052)
v1 = varMetadata(gse8052)
p2 = pData(gse8052_2)
v2 = varMetadata(gse8052_2)
pnew = rbind(p1,p2)
newp = new("AnnotatedDataFrame",data=pnew,varMetadata=v1)
gse8052@phenoData=newp
save(gse8052,file='GSE8052.eset')
write.table(pnew,row.names=F,sep='\t',file='GSE8052.pheno',quote=F)

"""

import sys, string, os, tempfile, subprocess, glob

rgeo="""loopGetGEO <- function(...) {
 returnGEO <- NULL
 notdone <- TRUE
 i = 0
 maxtry = 5
 while(notdone) {
   returnGEO <- try(getGEO(...))
   i = i+1
   if (i>maxtry) { 
      notdone <- FALSE
      if (inherits(returnGEO,"try-error")) {
         cat(paste('## loopGetGEO tried',i,'times - something bad going on - giving up'))
         return(returnGEO)
         } # bail
      } # i > maxtry
   if(!inherits(returnGEO,"try-error")) { notdone <- FALSE }
   else { cat("## loopGetGEO Caught an error in getGEO - maybe GEO ftp failure - so trying again....") }
   } # notdone
 return(returnGEO)
 } # loopGetGEO

loopGetGEOSuppFiles <- function(...) {
 returnGEO <- NULL
 notdone <- TRUE
 maxtry = 5
 i = 0
 while(notdone) {
  returnGEO <- try(getGEOSuppFiles(...))
  i = i+1
  if (i > maxtry) {
     notdone = FALSE
     if (inherits(returnGEO,"try-error")){
        cat(paste('## loopGetGEOSuppFiles tried',i,'times - something bad going on - giving up'))
        return(returnGEO)
       }
     } # i > maxtry
  if(!inherits(returnGEO,"try-error")) { notdone <- FALSE }
  else {cat("## loopGetGEOSuppFiles Caught an error in getGEOSuppFiles - maybe GEO ftp failure - so trying again....") }
  } # notdone
 return(returnGEO)
}



doLocal<-function(eset,outprefix,tdir)
{
    outf_pheno=paste(outprefix,'pheno',sep='.')
    outf_features=paste(outprefix,'features',sep='.')
    outf_eSet = paste(outprefix,'eset',sep='.') 
    phe = pData(eset)
    vm = varMetadata(eset)
    # need to get rid of :,;~| and other popular delimiters
    # these break , delimited parameter passing so must be dealt
    # whenever we save a file - no idea where you came from!
    rexphead = "[^(a-z,A-Z,0-9)]" # allowed in column values and names
    rexpval = "[^(a-z,A-Z,0-9,'.',':','/','_')]" # allowed in column values and names
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
    epath = file.path(tdir,outf_eSet)
    save(eset,file=epath)
    return(eset)
  }

rgeo=function(geoID,tdir)
{
# loopGetGEO is as suggested by stuart davis,
# the GEOquery BioC package maintainer
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
# Call as:  R CMD BATCH '--vanilla --slave --args GDS10 eSet maList ' geoquery.R logf
# both out paths will be supplied by Galaxy
# geoID will be supplied by user as a parameter
#
# note we need to convert GSE's differently from GDS geo datasets
# and for the moment, let's ignore that nasty complexity
res = ''
makeMAList = F
outf_eSet = file.path(tdir,paste(geoID,'eset',sep='.'))
outf_maList = file.path(tdir,paste(geoID,'malist',sep='.'))
cat('# geoquery.R, part of the Rexpression Galaxy toolkit http://esphealth.org/trac/rgalaxy')
cat(paste('# Starting at ',date()))
cat(paste('# geoquery.R got',geoID,'as the GEO id and',outf_eSet,'as outfile'))
library(GEOquery)
geo_prefix = substring(geoID,1,3) # eg GDS
print(paste('### geo prefix=',geo_prefix))
if (geo_prefix == 'GDS')
{
 # GEO DataSets (GDSxxx) are curated sets of GEO Sample data. A GDS record represents
 # a collection of biologically and statistically comparable GEO Samples and forms the 
 # basis of GEO's suite of data display and analysis tools. Samples within a GDS refer 
 # to the same Platform, that is, they share a common set of probe elements. 
 # Value measurements for each Sample within a GDS are assumed to be calculated 
 # in an equivalent manner, that is, considerations such as background processing 
 # and normalization are consistent across the dataset. Information reflecting 
 # experimental design is provided through GDS subsets.
 gds = loopGetGEO(geoID,AnnotGPL=T) 
 platrecgs = Meta(gds)$platform
 platrec = getGEO(platrecgs)
 # why bother with malist? deprecated I think
 if (makeMAList)
   {
     maList = GDS2MA(gds,GPL=platrec) # make an MAList (named maList when loaded!)
     save(maList,file=outf_maList) # supplied by Galaxy for the user's history
     print('maList saved')
   }
 eSet = GDS2eSet(gds) 
 # make an ExpressionSet (named eSet when loaded!)
 outprefix = geoID
 neweSet = doLocal(eSet,outprefix,tdir)
 print(paste(geoID,'eSet and pheno saved'))
 print(paste('# Done at ',date()))
} # GDS

if (geo_prefix == 'GSE')
{
# The GSE is the most confusing of the GEO entities. A GSE entry can represent an 
# arbitrary number of samples run on an arbitrary number of platforms. The GSE has 
# a metadata section, just like the other classes. However, it does not have a 
# GEODataTable. Instead, it contains two lists, accessible using GPLList and GSMList, 
# that are each lists of GPL and GSM objects. 
 i = 1
 esetList = loopGetGEO(geoID,GSEMatrix=T)
 dummy = loopGetGEOSuppFiles(geoID,makeDirectory=F,baseDir=tdir)
 if (length(esetList) > 1)  {
   comb = esetList[[1]]
   print(paste('### esetlist > 1, Annotation1 = ',comb@annotation))
   for (j in c(2:length(esetList)))
   {
     anotherEset = esetList[[j]]
     x = combine(comb,anotherEset)
     comb = x
     rm(x)
   }
   outprefix = geoID
   fixedEset = doLocal(eSet,outprefix,tdir)
   print(paste('### Saving, Annotation = ',fixedEset@annotation))
   print(paste(geoID,'Combined eSet and pheno saved'))
   print(paste('# Done at ',date()))
 
   i = 1
   for (e in esetList)
   {
    outf_eSet = file.path(tdir,paste(geoID,'_',i,'.eset',sep='')) # sequential in case > 1
    outprefix = paste(geoID,'_',i,sep='')
    fixedE = doLocal(e,outprefix,tdir)
    print(paste(geoID,'eSet and pheno saved'))
    print(paste('# Done at ',date()))
   }
  } 
  else {
    comb = esetList[[1]]
    outprefix = geoID
    fixedE = doLocal(comb,outprefix,tdir)
    print(paste('### Saving, Annotation = ',fixedE@annotation))
    print(paste(geoID,'eSet and pheno saved'))
    print(paste('# Done at ',date()))
 }
} # GSE

if (geo_prefix == 'GSM')
{
 # we could get the list of esets as esetlist = getGEO(gs,GSEMatrix=T)
 # but for now let's bail out
 print('### Sorry, this version does not handle GEO Series data sets')
} # GSM
if (geo_prefix == 'GPL')
{
 print('Sorry, this version does not handle GEO platform data sets')
} # GPL

}

# call our function
rgeo('%s','%s')
# instantiate with geoid and output directory
"""


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
    opath = os.environ.get('PATH', '').split(':')
    opath.insert(0,'/usr/local/bin/')
    opath.insert(1,'/usr/bin/')
    for path in opath:
        p = os.path.join(path, program)
        if (os.path.exists(p)) and (not os.path.isdir(p)):
            return p
    return None


def RRun(rcmd=[],outdir='.',title='myR',rexe=None):
    """
    run an r script, lines in rcmd,
    in a temporary directory
    move everything, r script and all back to outdir which will be an html file

    
      # test
      RRun(rcmd=['print("hello cruel world")','q()'],title='test')    
    """
    rlog = []
    flist = []
    if rexe == None:
        print '## bad news hello from rrun'
        rlog.append('## Error: Unable to find an executable R - sorry, please ask someone to fix this?')
        return rlog,flist # for html layout
    else:
        print '## rrun, rexe=%s' % rexe
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
    cruft = '.'*10 # from geoquery
    rlog = file(stofname,'r').readlines()
    rlog = [x.replace(cruft,'') for x in rlog]
    rlog = [x for x in rlog if x.strip() > '']# cruft
    flist = os.listdir(outdir)
    flist.sort
    flist = [(x,x) for x in flist] 
    for i,x in enumerate(flist):
        if x == rname:
            flist[i] = (x,'R script for %s' % title)
        elif x == stoname:
            flist[i] = (x,'R log for %s' % title)        
    return rlog,flist # for html layout
            
  

def main():
    """
    called as
    <command interpreter="python">
    reGeoquery.py '$gid'  '$title' '$htmlout' '$htmlout.extra_files_path' '${GALAXY_DATA_INDEX_DIR}/rg/bin/R'
    </command>
    """
    nparm = 4
    appname = os.path.split(sys.argv[0])[-1]
    if (len(sys.argv) < nparm):
        print '%s needs %d parameters - given %d as %s' % (appname,nparm,len(sys.argv),';'.join(sys.argv))
        sys.exit(1)
        # short command line error
    appname = sys.argv[0]
    geoid = sys.argv[1].strip()
    title = sys.argv[2].strip()
    logfname = sys.argv[3].strip()
    tdir = sys.argv[4].strip()
    try:
       os.makedirs(tdir)
    except:
       pass
    rexe = sys.argv[5]
    replace = string.whitespace + string.punctuation
    ttab = string.maketrans(replace,'_'*len(replace))
    title = string.translate(title,ttab)
    logf = file(logfname,'w')
    logf.write(galhtmlprefix % appname)
    logf.write('# %s - part of the Rexpression Galaxy toolkit http://rgenetics.org<br>\n' % (appname))
    logf.write('# Got parameters %s\n' % ' '.join(['"%s"' % x for x in sys.argv]))
    R = rgeo % (geoid,tdir)
    rlog,flist = RRun(rcmd=R,title=title,outdir=tdir,rexe=rexe) # run
    flist.sort()
    logf.write('<div>The following files were produced by this job<br>\n')
    logf.write('<table cellpadding="3" cellspacing="3">\n')
    for (fname,blurb) in flist:
        kind = '' # no type="foo"
        fpath = os.path.join(tdir,fname)
        size = ''
        if os.path.isfile(fpath):
            size = '(%1.1eMB)' % (float(os.path.getsize(fpath))/2**20)
        x = os.path.splitext(fname)[-1]
        if x in ['.log','.R']:
            kind = ' type="text/plain"'
        elif x == '.pheno':
            kind = ' type="application/vnd.ms-excel"'
        logf.write('<tr><td><a href="%s" %s>%s %s</a></td></tr>' % (fname,kind,blurb,size)) # or so
    logf.write('</table></div>\n')
    logf.write('<div>The R log from this run was: <br/><pre>\n')
    logf.write('\n'.join(rlog))
    logf.write('</pre></div>\n')
    logf.write(galhtmlpostfix)
    logf.write('\n')
    logf.close()

if __name__ == "__main__":
    main()
