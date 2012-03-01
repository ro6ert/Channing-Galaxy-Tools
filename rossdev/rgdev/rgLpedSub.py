#!/usr/bin/env python2.4

"""
Creates an html file to be used for viewing rgenetics files

    <command interpreter="python">
     rgLpedPbed.py '$lpedpath'  $log 
    </command>

Substantially updated October 2008 to allow multiple input lped files to be
selected from the default lped library and imported.

TODO add separate import of lped into history as well as the pbed datatype?


"""

import sys, os, glob, subprocess, shutil, tempfile

mng = '### makenewgalaxy' # used by exec_after_process hook to parse out new file paths/names

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
""" # this is used as the header for the rgenetics html file

def getMissval(inped=''):
   """ 
   read some lines...ugly hack - try to guess missing value
   should be N or 0 but might be . or -
   """
   commonmissvals = {'N':'N','0':'0','n':'n','9':'9','-':'-','.':'.'}
   try:
       f = file(inped,'r')
   except:
       return None # signal no in file
   missval = None
   while missval == None: # doggedly continue until we solve the mystery
        try:
          l = f.readline()
        except:
          break 
        ll = l.split()[6:] # ignore pedigree stuff
        for c in ll:
            if commonmissvals.get(c,None):
               missval = c 
               f.close()
               return missval
   if not missval: 
       missval = 'N' # punt
   close(f) 
   return missval

def doConvert(destdir=None,sourcepath=None,basename=None,flog=None,tlog=None,region='',rslist=''):
    """
echo 'Rgenetics http://rgenetics.org Galaxy Tools ped file importer'
echo "CL = $1 $2 $3 $4 $5"
echo 'called as importPed.sh $i $filelib_library_path $o $logfile $m $userId $userEmail'
echo "Calling Plink to import $1"
mkdir $4
echo 'Rgenetics http://rgenetics.org Galaxy Tools ped file importer'
echo "CL = $1 $2 $3 $4 $5"
echo 'called as importPed.sh $i $filelib_library_path $o $logfile $m $userId $userEmail'
echo "Calling Plink to import $1"
mkdir $4
cp $1/$2.map $4/
plink --file $1/$2 --make-bed --out $4/$2 --missing-genotype $5 > $3[rerla@hg rgenetics]$
    """
##    try:
##        os.makedirs(destdir)
##    except:
##        pass
    progname = os.path.split(sys.argv[0])[-1]
    inroot = os.path.join(sourcepath,basename)
    outroot = os.path.join(destdir,basename)
    mapsrc = '%s.%s' % (inroot,'map')
    mapdest = '%s.%s' % (outroot,'map')
    pedsrc = '%s.%s' % (inroot,'ped')
    peddest = '%s.%s' % (outroot,'ped')
    missval = getMissval(inped = pedsrc)
    if not missval:
        print >> sys.stdout, '### rgLpedPbed cannot read %s' % pedsrc
        missval = '0'
    if os.name == 'posix': # symlinks only available in posix systems
        os.symlink(mapsrc,mapdest)
        os.symlink(pedsrc,peddest)
    else: # otherwise copy - potentially expensive for large files
        shutil.copy(mapsrc,mapdest)
        shutil.copy(pedsrc,peddest)
    cl = ['plink --noweb --file',inroot,'--make-bed --out',outroot,
          '--missing-genotype',missval]
    #grab it all or add subset
    if region > '':
        try: # TODO make a regexp?
            c,rest = region.split(':')
            c = c.replace('chr','')
            rest = rest.replace(',','') # remove commas
            spos,epos = rest.split('-')
            spos = '%d' % int(float(spos)/1000.) # kb
            epos = '%d' % int(float(epos)/1000.) #kb
            s = '## %s parsing chrom %s from %s to %s\n' % (progname,c,spos,epos)
            flog.write(s)
            #plink --bfile mydata --chr 2 --from-kb 5000 --to-kb 10000
            cl += ['--chr',c,'--from-kb',spos,'--to-kb',epos]
        except:
            s = '##! %s unable to parse region %s - MUST look like "chr8:10,000-100,000' % (progname,region)
            print >> sys.stderr, s
            flog.write(s)
            flog.write('\n')
            flog.close()
            sys.exit(1)
    else:
        if rslist > '':
            cl += ['--snps',rslist.replace(' ','')]      
    cl  += ['> %s' % tlog]
    s = 'Calling plink with %s' % ' '.join(cl)
    flog.write(s)
    p = subprocess.Popen(' '.join(cl),shell=True,cwd=destdir)
    retval = p.wait() # run plink
    return missval
   

def doImport():
    """ convert lped into pbed and import into one of the new html composite data types for Rgenetics
        Dan Blankenberg with mods by Ross Lazarus 
        October 2007
        region or rslist added november 10 for ashg :)
        ross lazarus

         <command interpreter="python">
        rgLpedSub.py '$i' '$log' '$log.extra_files_path' '$region' '$rslist' 
        </command>

    """
    newbasenames = []
    progname = sys.argv[0]

    lp = sys.argv[1].replace('XX',',')
    logf = sys.argv[2].strip()
    tdir = sys.argv[3]
    region = sys.argv[4].strip() # gets precedence
    rslist = sys.argv[5].strip()
    try:
        os.makedirs(tdir)
    except:
        pass
    lpedpaths = list(lp.split(',')) # in case more than one
    flog = file(logf,'w')
    for p in lpedpaths: # for each file to be processed
        tfd,tlog = tempfile.mkstemp(prefix='rgLpedSublog')
        importpath,basename = os.path.split(p)
        repository,rest = os.path.split(importpath) # get repository root
        missval = doConvert(destdir=tdir,sourcepath=importpath,basename=basename,flog=flog,
            tlog=tlog,region=region,rslist=rslist)
        plinklog = file(tlog,'r').read()
        os.unlink(tlog)
        flog.write('Converted %s, using missing value code = %s\n' % (p,missval))
        flog.write(plinklog)
        flog.write('\n')
        outfile = os.path.join(tdir,'%s.html' % basename)
        out = file(outfile,'w') # start again
        out.write(galhtmlprefix % progname)
        out.write('<h2>Rgenetics data imported from %s</h2>\n<hr><ul>' % basename)
        flist = glob.glob(os.path.join(tdir, '%s.*' % (basename))) # all the files matching basename
        for i, data in enumerate( flist ):
            ofn = os.path.split(data)[-1]
            out.write('<li><a href="%s">%s</a></li>\n' % (ofn,ofn))
        out.write('<hr><h3>Plink log follows:</h3><pre>')
        out.write(plinklog)
        out.write('</pre><br>')
        out.write("</div></body></html>")
        out.close()
        newbasenames.append(basename)
    for b in newbasenames:
        s = '%s\t%s\t%s\n' % (mng,tdir,b)
        flog.write(s) # used by exec_after_process to make new datasets
        print >> sys.stdout, s
    flog.close()


if __name__ == "__main__": 
   doImport()

