""" 
reUpdate_eset.py
Sometimes, the phenodata supplied have characters that confuse our
dumb parameter passing mechanisms - eg , and sometimes BioC things
use : or ; as delimiters.
We clean these up when we bring data in but this utility is provided in case
you have problems

copyright ross lazarus
sept 2009
licensed under the terms of the LGPL
for the rexpression suite
    
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
#( infile,outdir,eoutpath,poutpath)
rprog="""sanitiseChars=function(infile="%s",tdir="%s",eoutpath="%s",poutpath="%s",
        newppath="%s", newanno="%s")
{
    library('Biobase')
    dname = load(infile)
    eset = get(dname) # R magic - we don't actually know the objects' R identity...
    if (newppath != 'None') # replace phenodata
    {
        phe = read.table(file=newppath,sep='\t',header=T)
        vm <- data.frame(labelDescription = colnames(phe)) # update varMetadata slot
    }
    else 
    {
        phe = pData(eset)
        vm = varMetadata(eset)
    }

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
    if (newanno > "")
    {
      eset1@annotation=newanno
    }
    save(eset,file=eoutpath)
    write.table(phe,file=poutpath,sep='\t',row.names=F,quote=F)
}
# call
sanitiseChars()"""

## R = sanitiseChars % ( infile,outdir,eoutpath,poutpath,newpheno,newanno)


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


def RRun(rcmd=[],outdir='.',title='myR', rexe=None):
    """
    run an r script, lines in rcmd,
    in a temporary directory
    move everything, r script and all back to outdir which will be an html file

    
      # test
      RRun(rcmd=['print("hello cruel world")','q()'],title='test')    
    """
    rlog = []
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
reUpdate_eset.py "$i.extra_files_path" "$title" "$outhtml" "$outhtml.files_path"
 "$i.metadata.base_name" "$i.ext" "R" "$p" "$c" "$anno"
</command>
    writes a new rexpression object filtered by phenovals present in phenocol
    """
    nparm = 9
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
    outext = sys.argv[6].strip()    
    rexe = sys.argv[7].strip()
    newpheno = sys.argv[8].strip()
    newpid = sys.argv[9].strip() # need the id col if not the first
    newanno = sys.argv[10].strip() # need the id col if not the first
    replace = string.whitespace + string.punctuation
    ttab = string.maketrans(replace,'_'*len(replace))
    title = string.translate(title,ttab)
    newhtml = [galhtmlprefix % title,] 
    newhtml.append('<div class="document"># %s - part of the Rexpression Galaxy toolkit http://esphealth.org/trac/rgalaxy\n' % (appname))
    newhtml.append('# Got parameters %s\n' % ' '.join(sys.argv))
    newhtml.append('Files created by this job:')
    newhtml.append('<ol>\n')
    if newpheno <> 'None': # check that id's are compatible
       try:
          newpid = int(newpid)
       except:
          newhtml.append('<h3>Error - new pheno file id column = %s</h3>' % newpid)
          newpid = 1 # punt
       npd = file(newpheno,'r').readlines()
       npd = [x.strip().split('\t') for x in npd]
       oldpheno = os.path.join(infpath,'%s.pheno' % basename)
       try:
          opd = file(oldpheno,'r').readlines()
          opd = [x.strip().split('\t') for x in opd]
       except:
          print >> sys.stderr, '## sorry, unable to open old pheno file %s' % (oldpheno)
          print >> sys.stderr, '## cannot complete job %s' % title
          sys.exit(1)     
       if newpid >= 1:
           c = newpid - 1
           newids = [x[c] for x in npd]
       else:
	   newids = []
       mismatch = ['%s <> %s' % (x,opd[i][0]) for i,x in enumerate(newids) \
           if (i > 0) and (x <> opd[i][0])] # ignore first row
       if len(mismatch) > 0:
          newhtml.append('<div><h3>WARNING - found mismatching ids</h3>')
          newhtml.append('<table>')
          for i,m in enumerate(mismatch):
              newhtml.append('<tr><td>%d</td><td>%s</td></tr>' % (i,m))
          newhtml.append('</table></div>')
       if newpid > 1: # need to rewrite with id as first column
          c = newpid -1
          dummy, tpheno = tempfile.mkstemp(prefix=title)
          npd = [[x[c],] + x[:c] + x[(c+1):] for x in npd] # move id to first column
          npd = ['\t'.join(x) for x in npd]
          f = file(tpheno,'w')
          f.write('\n'.join(npd))
          f.write('\n')
          f.close()
          newpheno = tpheno # so R will use this
    # and convert into R code to make a vector
    try:
        os.makedirs(outdir)
    except:
        pass # already exists
    exts = ['.affybatch','.eset','.malist']
    inflist = glob.glob(os.path.join(infpath,'*.*'))
    validinfiles = [x for x in inflist if os.path.splitext(x)[-1] in exts]
    infile = None
    if len(validinfiles) > 0:
        infile = validinfiles[0]
    if not infile:
        print '%s requires an input file - supplied one, no %s.eset/.affybatch/.malist could be opened' % (appname,infpath)
        sys.exit(1) 
    eoutpath = os.path.join(outdir,'%s.eset' % basename)
    poutpath = os.path.join(outdir,'%s.pheno' % basename)  
    R = rprog % (infile,outdir,eoutpath,poutpath,newpheno,newanno)
    # must restore real phenotype codes - encoded to hide commas in command line pass
    flist,rlog = RRun(rcmd=R,outdir=outdir,title=title,rexe=rexe)
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
    try:
        os.unlink(tpheno)
    except:
        pass

if __name__ == "__main__":
    main()




