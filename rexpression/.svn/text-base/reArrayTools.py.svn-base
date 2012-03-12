"""
we need to make html links to genbank and such from the xls we
generate if there's no db annotation available
Unfortunately, then we have no idea what columns to expect in the featureData slot
since it's a zoo out there in arrayexpress and geo land. 

reArrayTools.py from
reAffyExpress.py
A Galaxy tool based on the simplest dichotomous example in AffyExpress
documentation I could find. Use this as a prototype to write more
complex models.
This code called by the job runner, parameters as shown in the corresponding
xml tool descriptor.

From ?regress

regress             package:AffyExpress             R Documentation

Run regression to fit genewise linear model

Description:

     Fit genewise linear model using LIMMA package, ordinary linear
     regression,  or permutation method.

Usage:

     regress(object, design, contrast, method, adj="none", permute.time=1000) 

Arguments:

  object: an "ExpressionSet"

  design: design matrix from the make.design function

contrast: contrast matrix from the make.contrast function

  method: Three methods are supported by this function: "L" for using
          LIMMA method - compute moderated t-statistics and log-odds 
          of differential expression by empirical Bayes shrinkage of
          the standard  errors towards a common value,   "F" for using
          ordinary linear regression,  "P" for permuation test by
          resampling the phenotype

     adj: adjustment method for multiple comparison test, including
          "holm",  "hochberg", "hommel", "bonferroni", "BH", "BY",
          "fdr", "none".  The default value is "none". Type
          help(p.adjust) for more detail.

permute.time: number of permutation times, only used for the
          permutation  method.

:
Value:

     A dataframe contains rows for all the genes from object and the
     following columns: ID(probeid); Log2Ratio (estimate of the effect
     or the contrast,  on the log2 scale);  F (F statistics); P.Value
     (raw p-value); adj.P.Value (adjusted  p-value or q-value)

Author(s):

     Xiwei Wu xwu@coh.org, Xuejun Arthur Li xueli@coh.org

References:

     Smyth, G.K. (2005) Limma: linear models for microarray data. In:
     Bioinformatics and Computational Biology Solutions using R and 
     Bioconductor, R. Gentleman,  V. Carey, S. Dudoit, R. Irizarry,  W.
     Huber (eds.), Springer, New York, pages 397-420

Examples:

     data(testData)
     normaldata<-pre.process("rma",testData)

     ## Create design matrix
     design<-make.design(pData(normaldata), "group")

     ## Create contrast matrix - Compare group "A" vs. "C"
     contrast<-make.contrast(design, "A", "C")

     ## Identify differentially expressed gene by using LIMMA method
     result<-regress(normaldata, design, contrast, "L")


 AffyQA(parameters = c("estrogen", "time.h"), raw = raw)

  result.wrapper <- AffyRegress(normal.data = filtered, cov = "estrogen",
+     compare1 = "present", compare2 = "absent", method = "L",
+     adj = "fdr", p.value = 0.05, m.value = log2(1.5))

 target<-data.frame(drug=(c(rep("A",4),rep("B",4),rep("C",4))), 
     gender=factor(c(rep("M",6),rep("F",6))), 
     group=factor(rep(c(1,2,3),4)))

     # Example1:  Compare drug "A" vs. "B"
     design1<-make.design(target, "drug")
     contrast1<-make.contrast(design1, "A", "B")

     # Example2:  Compare drug "A" vs. "B", adjusting for "group" variable
     design2<-make.design(target, c("drug","group"))
     contrast2<-make.contrast(design2, "A", "B")

     # Example3:  Suppose you are interested in "drug", "group" interaction
     design3<-make.design(target, c("drug","group"), int=c(1,2))
     contrast3<-make.contrast(design3, interaction=TRUE)

     # Example4:  Compare drug "A" vs. "B" among "male"
     # Notice that you must use an design matrix containing an interaction term
     design4<-make.design(target, c("drug","gender"), int=c(1,2))
     contrast4<-make.contrast(design4, "A", "B", "M")

    
There are serious problems with annotation - affy seems to work but illumina is less tractable
From bioc,
Dear Seth,
> Thanks for the useful suggestions. I'm digging the source code of the
> annotate package, it seems really interesting...
> Below a chunk of my code that use require:
>
> imp.data<- function(...){
>       # importing raw.data
>       ...
>       ...
>       ## load the proper annotation package regardless the technology
>       if( !grepl("\\.db", raw.data@annotation) ) {
>           raw.data@annotation<- paste(raw.data@annotation, ".db", sep="")
>       }
>       require( raw.data@annotation, character.only=T )
>
>       return(raw.data)
> }
>
> raw.data is an object of class 'ExpressionSet'.
>
> Thanks again for the help!
>
> Regards,
> Paolo
 
"""
import  sys, string, os, tempfile, shutil, time, subprocess
 

oursep = '=+comma+=' # shouldn't be in any data?

galhtmlprefix = """<?xml version="1.0" encoding="utf-8" ?>
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
<head>
<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
<meta name="generator" content="Galaxy tool output - see http://galaxyproject.org/" />
<title>Tool is %s</title>
<link rel="stylesheet" href="/static/style/base.css" type="text/css" />
</head>
<body>
<div class="document">
"""

galhtmlpostfix = """</div></body></html>"""
# eg
# library(AffyExpress)
# load( "Golub_Merge.eset")

#  reFunc('test.log', "Golub_Merge.eset" , "golub", '.','eset','rma','Gender','F','M','L','fdr','0.001','1.5','test.html',0)
#  reFunc('test.log', "ETABM25.affyBatch" , "etabm25", '.','eset','rma','FactorValue..Age.','<40','>40','L','fdr','0.001','1.5','test.html',0)

##> slotNames(select)                                                                      
## [1] "ID"                      "foldChange"                                              
## [3] "FValue"                  "pValue"                                                  
## [5] "adjPVal"                 "contrast"                                                
## [7] "regressionMethod"        "adjustment"                                              
## [9] "significantIndex"        "significantPvalueCutoff"                                 
##[11] "significantFCCutoff"     "fileName"                                                
##[13] "annotation"              "normalizationMethod"                                     
##[15] "filterMethod"            ## fd = eSet@featureData                                                              
##> slotNames(fd)                                                                     
##[1] "varMetadata"       "data"              "dimLabels"                             
##[4] ".__classVersion__"                                                             
##> fdd = fd@data                                                                     
##> names(fdd)                                                                        
## [1] "ID"                    "Gene.title"            "Gene.symbol"                  
## [4] "Gene.ID"               "UniGene.title"         "UniGene.symbol"               
## [7] "UniGene.ID"            "Nucleotide.Title"      "GI"                           
##[10] "GenBank.Accession"     "Platform_CLONEID"      "Platform_ORF"                 
##[13] "Platform_SPOTID"       "Chromosome.location"   "Chromosome.annotation"        
##[16] "GO.Function"           "GO.Process"            "GO.Component"                 
##[19] "GO.Function.1"         "GO.Process.1"          "GO.Component.1"               
##>                            

reFunc="""
#Function to compute a design matrix for making all possible pairwise comparisons

design.pairs <- function(levels) {
n <- length(levels)
design <- matrix(0,n,choose(n,2))
rownames(design) <- levels
colnames(design) <- 1:choose(n,2)
k <- 0
for (i in 1:(n-1))
for (j in (i+1):n) {
k <- k+1
design[i,k] <- 1
design[j,k] <- -1
colnames(design)[k] <- paste(levels[i],"-",levels[j],sep="")
}
design
}




reFunc <- function(infpath='%s',infname='%s',tDir='%s',
inftype='%s',norm_method='%s',phename='%s',pheval1='%s',
pheval2='%s',meth='%s',padj='%s',pthresh=%s,mthresh=%s,topTablename='%s',
permutations=%s)
{ 
library("ArrayTools")
library("affy")
sessionInfo()
x = load(infpath)
alleset = get(x)
alleset = na.omit(alleset)
klass=class(alleset)
eset = alleset[, ((alleset[[phename]]==pheval1 | alleset[[phename]]==pheval2) & !is.na(alleset[[phename]]))]
# remove all other values or arraytools design/contrast matrix code fails.
# TODO fixme!!
if (inftype == 'affybatch')
    {
    # this makes AffyQA.html and some low res png's 
    # TODO MAKE PDF's or better resolution so the labels are legible for complex data
    AffyQA(parameters = phename, raw = eset)
    normaldata <- pre.process(method = norm_method, raw = eset)
    eset <- normaldata
    eset = na.omit(eset)
    }
## Create design matrix
# remove all NA's or things fail
#print('Removing NAs from eset')
#eset<-eset[,is.na(exprs(eset)[1,])==F]
#eset<-eset[is.na(exprs(eset)[,1])==F,]
print('Getting pheno data from eset')
phe = pData(eset)
print(paste('Getting phename vector',phename,'from pheno data'))
p = phe[,phename]
up = unique(p)
print(paste('Unique values in',phename,'are',up))
if (length(up) > 20) { 
   # arbitrary
   print(paste('### ',length(up),' is too many levels to be worth computing IMHO'))
  }
print(paste('Phenotype chosen =',phename,' '))                                            
design<-new("designMatrix",target=phe,covariates=phename)                                
print(paste('design matrix for',phename,'='))                                            
print(design)                                                                                   
contrast<-new("contrastMatrix",design.matrix=design,compare1=pheval1,compare2=pheval2)   
print(paste('contrast for',phename,':',pheval1,'vs',pheval2,'='))                        
print(contrast)
if (meth == 'permutation') {                                                               
    result<-regress(eset, contrast=contrast, method=meth, adj=padj, permute.time=permutations) }
    else {
    result<-regress(eset, contrast=contrast, method=meth, adj=padj) }
print(summary(result))
select <- selectSigGene(result, p.value = pthresh, fc.value = log2(mthresh))              
                                                                  
mydir=getwd()

if (sum(getIndex(select),na.rm=T) == 0) {                                                               
    print("####No significant results were found. You may wish to try again with a larger P value threshold?") }
else { 
    # worth trying to produce outputs                                                                               
    gotAnno = F
    a = getAnnotation(select)
    if (length(a) > 0) {
     gotAnno = T
     packageName <- paste(a, ".db", sep="")
     if (! require(packageName,character.only=T))
        {
           print(paste('### annotation package',packageName,'must be installed to proceed'))
           source('http://bioconductor.org/biocLite.R')
           biocLite(packageName)
        }
        if (! require(packageName,character.only=T)) {
           print(paste('### annotation package',packageName,'could not be auto-installed'))
           print(paste('### please get someone to install it manually?'))
           gotAnno = F }

     } # length a

    if ((klass != 'LumiBatch') & gotAnno) {
         # arraytools annotation seems to barf on NBTcamp299 lumibatch  
        library(packageName , character.only=TRUE)                                            
        require("annaffy")                              
        probeset_id  <- getID(select)[getIndex(select)]
        print(paste('## using',packageName,'probeset_id[1:10]=',probeset_id[1:10]))
        selected = getIndex(select)
        FC <- as.data.frame(getFC(select))                                                    
        FC.select <- FC[getIndex(select),]                                            
        P.Value <- getP(select)[getIndex(select)]                                             
        Adj.P.Value <- getAdjP(select)[getIndex(select)]                                      
                                                                                              
        if (sum(is.na(P.Value)) == length(P.Value)){                                          
            fullresult <- data.frame(Probeset_id  = probeset_id , FC) } 
            else {                                                                          
            fullresult <- data.frame(Probeset_id  = probeset_id , FC.select,P.Value =  P.Value, Adj.P.Value = Adj.P.Value)}
        anncols <- aaf.handler()[c(2,3,4,6,7,8,9,10)]                                         
        anntable <- aafTableAnn(probeset_id, packageName , anncols)                           
        testtable <- aafTable(items=fullresult)
        otable = merge(testtable,anntable)
        sortedtable = otable[order(Adj.P.Value),]
        saveText(sortedtable,file=paste(topTablename,'xls',sep='.'))
        saveHTML(sortedtable, paste(topTablename, "html", sep = "."), title = paste("Differentially Expressed Genes by",phename))
        
    }
    else { 
    # try featureData simple output - no html fancy stuff
    fd = eset@featureData                                                                 
    fdd = fd@data                                                                         
    if (length(fdd$ID) > 0) {
        probeset_id  <- getID(select)[1:length(getP(select))]                                 
        FC <- as.data.frame(getFC(select))                                                    
        selected <- getIndex(select)                                                          
        P.Value <- getP(select)                                                               
        Adj.P.Value <- getAdjP(select)
        #Genb = fdd[,10]
        # use it all: was anno = fdd[,c(1,2,3,4,5,6,7,9)]                                                      
        if (sum(is.na(P.Value)) == length(P.Value)) {                                          
            fullresult <- data.frame(GenBankAcc=Genb, Probeset_id  = probeset_id ,FoldChange=FC) }
        else {                                                                          
            fullresult <- data.frame(ID=fdd[,1],
            logAdjP = -log10(Adj.P.Value),                             
            log10P =  -log10(P.Value) ,
            Probeset_id  = probeset_id ,FoldChange=FC,selected = selected) }
        allres=cbind(fullresult,fdd)                                                     
        res = subset(allres,selected==T)  
        rres = res[order(-res$logAdjP,res$ID),]                               
        write.table(rres,file=paste(topTablename,'xls',sep='.'),row.names=F,sep='\t',quote=F) }
    else  {
        print('### Unable to find an annotation source in your expression data eset - sorry') }
    } # simple output
  } # no sig results
} # end 
reFunc() # will run with default arguments
"""
## prepare script
## R = reFunc % (infpath,infname,outpath,inftype,'rma', phename,phevals[0].strip(),phevals[1].strip(),
## meth,padj,pthresh,mthresh,topTablename,permutations)

debug = False # useful for debugging RRun()

def timenow():
    """return current time as a string
    """
    return time.strftime('%d/%m/%Y %H:%M:%S', time.localtime(time.time()))

def whereis(program):
    opath = os.environ.get('PATH', '').split(':')
    opath.insert(0,'/usr/local/bin/')
    opath.insert(0,'/usr/bin/')
    for path in opath:
        p = os.path.join(path, program)
        if (os.path.exists(p)) and (not os.path.isdir(p)):
            return p
    return None


def RRun(rcmd=[],outdir='.',title='myR', rexe=None):
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
    else:
        f.write(rcmd) # assume is string already
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
            flist[i] = (x,'R script for %s' % title)
        elif x == stoname:
            flist[i] = (x,'R log for %s' % title)        
    return rlog,flist # for html layout
            
 

def main():
    """called as

<command interpreter="python">
reArrayTools.py '$i.extra_files_path' '$title' '$outhtml' '$i.metadata.base_name'
 '$i.extension' '$phename' '$method' '$permutations' '$padj' '$pthresh' '$mthresh'
 '$outhtml.files_path' "$phevals" '${GALAXY_DATA_INDEX_DIR}/rg/bin/R'
</command>

    """
    sep = ','
    nparm = 14
    appname = sys.argv[0]
    if (len(sys.argv) < nparm):
        print '%s needs %d parameters - given %d as %s' % (appname,nparm,len(sys.argv),';'.join(sys.argv))
        sys.exit(1)
        # short command line error
    appname = sys.argv[0]
    infp = sys.argv[1].strip()
    title = sys.argv[2].strip()
    outhtml = sys.argv[3].strip() 
    infbasename = sys.argv[4].strip()
    inftype = sys.argv[5].strip()
    infpath = os.path.join(infp,'%s.%s' % (infbasename,inftype)) # build extension
    phevals = sys.argv[13].strip().split(',') # columns for design/contrast
    phename = sys.argv[6].strip() # first is name
    meth = sys.argv[7].strip() # L,F or P
    permutations = sys.argv[8].strip()
    try:
        p = int(permutations)
    except:
        permutations = '0'
    padj = sys.argv[9].strip() # fdr etc
    pthresh = sys.argv[10].strip()
    mthresh = sys.argv[11].strip()
    try:
        p = float(pthresh)
    except: 
        pthresh = '1.0' # all
    try:
        m = float(mthresh)
    except:
        mthresh = '0.0' # all
    outpath = sys.argv[12].strip() 
    try:
        os.makedirs(outpath)
    except:
        pass
    rexe = sys.argv[14].strip() 
    norm_method='rma'    
    replace = string.whitespace + string.punctuation
    ttab = string.maketrans(replace,'_'*len(replace))
    title = string.translate(title,ttab)
    ourlog = []
    ourlog.append('# %s - part of the Rexpression Galaxy toolkit http://esphealth.org/trac/rgalaxy\n' % (appname))
    ourlog.append('# Job %s Got parameters %s\n' % (title,' '.join(sys.argv)))
    if len(phevals) <> 2: #
        ourlog.append('### WARNING - phevals must have exactly 2 phenotype values to contrast - supplied with %d in %s' % (len(phevals),phevals) ) 
    try:
        os.makedirs(outpath) # not made yet?
    except:
        pass
    html = []
    html.append(galhtmlprefix % appname)
    html.append('<h2>Output from job title %s run at %s</h2><ol>' % (title,timenow()))
    topTablename = '%s_%s_Anova' % (phename,infbasename)
    s = ' '.join(('Top table is for',infbasename,'phenotype=',phename,'contrast=',phevals[0],
                  'vs',phevals[1],'method=',meth,'pthresh=',pthresh,'foldthresh=',mthresh,
                  'adjustment=',padj))
    html.append('topTable is for %s<br/>\n' % s)
    R = reFunc % (infpath,infbasename,outpath,inftype,'rma',
        phename,phevals[0].strip().replace(oursep,','),phevals[1].strip().replace(oursep,','),meth,padj,pthresh,mthresh,
        topTablename,permutations)
    ## TODO: fixme - we use a crude comma command line quoting system since otherwise we die if there are any
    ## in phenotype values - 
    rlog,flist = RRun(rcmd = R,outdir = outpath,title=topTablename,rexe=rexe)
    ourlog.append('## R log:\n')
    ourlog += rlog
    fl = os.listdir(outpath)
    gothtml = False
    outxls = False
    ourlog.append('first fl = %s' % fl)
    for fn in fl:
        fext = os.path.splitext(fn)[-1]
        if fext == '.html':
            gothtml = fn
        elif fext == '.xls':
            outxls = os.path.join(outpath,fn) # full path
    if not gothtml: # must rewrite our xls file
        if outxls:
            ourlog.append('## writing an html equivalent of %s' % outxls)
            outh = os.path.splitext(outxls)[0]
            outh = '%s.html' % outh
            outhf = file(outh,'w')
            fl = file(outxls,'r').readlines()
            fl = [x.split('\t')[:-1] for x in fl] # drop last \n without removing trailing tabs
            head = fl[0]
            try:
                i = head.index('selected')
            except:
                i = None
            if i:
                for n,row in enumerate(fl):
                    rowc = row[:(i-1)] + row[i:] # eliminate redundant column
                    fl[n] = rowc
            h = [galhtmlprefix % appname,]
            h.append('<table cellpadding="1" cellspacing="0" border="1">')
            for i,row in enumerate(fl):
                for j,x in enumerate(row):
                    if len(x) == 0:
                        row[j]="&nbsp;"
                        # trick so all cell borders appear...
                if i == 0: # head
                    s = "</th><th>".join(row)
                    h.append('<tr><th>%s</th></tr>\n' % s)
                else:
                    s = "</td><td>".join(row)
                    #s = "</td><td align='center' valign='center'>".join(row)
                    h.append('<tr><td>%s</td></tr>\n' % s)
            h.append('</table>')
            outhf.write('\n'.join(h))
            outhf.write('\n')
            outhf.close()
        else:
            ourlog.insert(0,'## WARNING: No html or xls file found - probably an error - sorry - clues may follow')
    ourlog.append('R code follows:\n%s\n' % R)
    fl = os.listdir(outpath)
    fl.sort()
    for f in fl:
        html.append('<li><a href="%s">%s</a></li>' % (f,f))
    html.append('</ol>')
    html.append('<hr><h2>Tool was %s. Job log Follows:</h2><hr/><br/><pre>%s</pre>' % \
                (appname,'\n'.join(ourlog))) # show log    html,write('\n'.join(ourlog))
    html.append(galhtmlpostfix)
    outf = file(outhtml,'w')
    outf.write('\n'.join(html))
    outf.write('\n')
    outf.close()

if __name__ == "__main__":
    main()

