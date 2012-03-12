"""
reArraySubset.py
Simple way of subsetting an expression object by some arbitrary phenotype column values
eg to keep some specific chips or phenotypes..
copyright ross lazarus
sept 2009
licensed under the terms of the LGPL
for the rexpression suite

some recipes in
http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/R_BioCondManual.html#biocon_affypack


# Filtering by expression values

      get.array.subset(x.rma, "treatment",c("CMP1","Ctrl")) # When R loaded *.CEL files it also reads experiment layout from covdesc.txt file (Ctrl and CMP1). The "get.array.subset" function makes it easy to select subsets of arrays.
      results <- pairwise.comparison(x.rma, "treatment", c("1", "2"), raw.data) # computes mean intensities of replicates from two treatments (means), calculates log2-fold changes between means (fc), performs t-test (tt) and writes out PMA calls (calls). Type "?pairwise.comparison" for more help on this function.
      write.table(data.frame(means(results), fc(results), tt(results), calls(results)), file="my_comp.txt", sep="\t") # Command to export all 'pairwise.comparison' data into one table.
      sort(abs(fc(results)),decreasing=TRUE)[1:100] # Prints out 100 highest fold-changes (fc), for means of intensities (means), for ttest (tt), for PMA (calls).
      significant <- pairwise.filter(results, min.exp=log2(10), min.exp.no=2, fc=log2(8), min.present.no=4, tt= 0.001, present.by.group=FALSE) # Function 'pairwise.filter' takes the output from 'pairwise.comparison' and filters for significant changes. Filter Arguments:
          o min.exp: minimum expression cut off
          o min.exp.no: occurence of 'min.exp' in at least this number of chips
          o min.present.no: present calls on at least this number of chips
          o present.by.group: If true, then present count is restricted to replicate groups and not all chips of an experiment!
          o fc: A gene must show a log2 fold change greater than this to be called significant. A 2.5-fold change can be specified like this: 'fc=log2(2.5)'
          o tt: A gene must be changing with a p-score less than this to be called significant
      Type "?pairwise.filter" for more info on the different arguments of this function.

      write.table(data.frame(means(significant), fc(significant), tt(significant), calls(significant)), file="my_comp.txt", col.names = NA, sep="\t") # Exports all 'pairwise.filter' data into one table.

# Plotting results

      plot(significant, type="scatter") # Plots significant changes from above as scatter plot. Alternative plotting types: type="ma" or type="volcano".
      plot(results,significant, type="scatter") # Plots means of the two replicate groups as scatter plot and high-lights all genes that meet criteria of 'significant' filter. Meaning of colors: red - all present, orange - all present in one group or the other, yellow - all that remain.
      plot(results,type="scatter") # Plots means of the two replicate groups as scatter plot.
      png(file="figure1.png"); plot(significant,type="scatter"); dev.off() # Writes scatter plot to image file.

on a raw (affybatch eg) dataset
here as a string called rprog

Challenge is to present the phenotype data in the xml so columns can be
specified.     
    
"""
import  sys, string, os, tempfile, subprocess, glob

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
oursep = '~*comma!~' # shouldn't be in any data?

rprog="""arraysubset=function(infile,phename,phevals,tdir,eoutpath,poutpath)
{
    library('simpleaffy')
    dname = load(infile)
    eset1 = get(dname) # R magic - we don't actually know the objects' R identity...
    subset1 = get.array.subset(eset1, phename,phevals)
    save(subset1,file=eoutpath)
    p = pData(subset1)
    write.table(p,file=poutpath,sep='\t',row.names=F,quote=F)
}
# call
arraysubset("%s","%s",c(%s),"%s","%s","%s")"""

## R = rprog % (infpath,phename,phevals,outdir,eoutpath,poutpath)


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


def RRun(rcmd=[],outdir='.',title='myR',rexe=None):
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
    else:
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
    flist.sort()
    for i,x in enumerate(flist):
        if x[1] == rname:
            flist[i] = (x,'## R script for %s' % title)
        elif x[1] == stoname:
            flist[i] = (x,'## R log for %s' % title)   
        else:
            flist[i] = (x,x)
    return flist,rlog # for html layout
            

def main():
    """called as 
<command interpreter="python">
reArraySubset.py '$i.extra_files_path' '$title' '$outhtml' '$outhtml.files_path'
 '$i.metadata.base_name' '$phename' '$phevals' '${GALAXY_DATA_INDEX_DIR}/rg/bin/R'
</command>
    writes a new rexpression object filtered by phenovals present in phenocol
    """
    nparm = 8
    appname = sys.argv[0]
    if (len(sys.argv) < nparm):
        print '%s needs %d parameters - given %d as %s' % (appname,nparm,len(sys.argv),';'.join(sys.argv))
        sys.exit(1)
        # short command line error
    appname = sys.argv[0]
    infpath = sys.argv[1].strip()
    title = sys.argv[2].strip()
    logfname = sys.argv[3].strip()
    outdir = sys.argv[4].strip()
    basename = sys.argv[5].strip()
    phename = sys.argv[6].strip() # allow multiples
    phevals = sys.argv[7].strip().split(',')
    rexe = sys.argv[8].strip()
    phevals = ','.join(["'%s'" % x for x in phevals]) # quote
    # and convert into R code to make a vector
    replace = string.whitespace + string.punctuation
    ttab = string.maketrans(replace,'_'*len(replace))
    title = string.translate(title,ttab)
    try:
        os.makedirs(outdir)
    except:
        pass # already exists
    infpath = os.path.join(infpath,basename)
    infile = '%s.affybatch' % infpath
    if not os.path.exists(infile):
        infile = '%s.eset' % infpath
        if not os.path.exists(infile):
            infile = '%s.malist' % infpath
            if not os.path.exists(infile):
                print '%s requires an input file - supplied one, %s does not seem to exist' % (appname,infpath)
                sys.exit(1) 
    eoutpath = os.path.join(outdir,'%s_subset.eset' % basename)
    poutpath = os.path.join(outdir,'%s_subset.pheno' % basename)    
    print '\n'.join((infile,phename,phevals,outdir,eoutpath,poutpath))
    R = rprog % (infile,phename,phevals.replace(oursep,','),outdir,eoutpath,poutpath)
    # must restore real phenotype codes - encoded to hide commas in command line pass
    flist,rlog = RRun(rcmd=R,outdir=outdir,title=title,rexe=rexe)
    newhtml = [galhtmlprefix % title,] 
    newhtml.append('<div class="document"># %s - part of the Rexpression Galaxy toolkit http://esphealth.org/trac/rgalaxy\n' % (appname))
    newhtml.append('# Got parameters %s\n' % ' '.join(sys.argv))
    newhtml.append('Files created by this job:')
    newhtml.append('<ol>\n')
    for f in flist:
        newhtml.append('<li><a href="%s">%s</a></li>' % (f[0],f[1]))
    newhtml.append('</ol></br></div>\n')
    newhtml.append('<div class="document">R log follows:<br/><pre>\n')
    newhtml += rlog # .append(''.join(rlog))
    newhtml.append('</pre></div>\n')
    newhtml.append(galhtmlpostfix)
    logf = file(logfname,'w')
    logf.write('\n'.join(newhtml))
    logf.write('\n')
    logf.close()


if __name__ == "__main__":
    main()

