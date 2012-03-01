# read linkage ped to snpmatrix
# ross lazarus may 2008 for the rgenetics project
# copyright ross lazarus 2008 
# released under the LGPL v2 or later.
"""
grrrr - if I call library('snpMatrix'), R writes 
Loading required package: survival
Loading required package: splines
Loading required package: hexbin
Loading required package: grid
Loading required package: lattice
to stderr - not stdout!!!

So Galaxy assumes the tool returned an error!

Took me a day to figure this out...

read.snps.pedfile Read genotype data from a LINKAGE "pedfile"

Description
This function reads data arranged as a LINKAGE "pedfile" with some restrictions and returns a list
of three objects: a data frame containing the initial 6 fields giving pedigree structure, sex and disease
status, a vector or a data frame containing snp assignment and possibly other snp infomation, and
an object of class "snp.matrix" or "X.snp.matrix" containing the genotype data

Usage
read.snps.pedfile(file, snp.names=NULL, assign=NULL, missing=NULL, X=FALSE, sep=".",..)
"""
import sys,os,string,tempfile,subprocess

progname = os.path.split(sys.argv[0])[-1]

library = '/usr/local/galaxy/data/rg'

galhtmlprefix = """<?xml version="1.0" encoding="utf-8" ?>
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" 
"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
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
# save a copy of the pedfile as a galaxy snpmatrix file
# grr - needs haploview style info rather than plink style map
# easy to fix - should take a .map file if available!

def doImport(outfile='test',flist=[],log=[]):
    """ import into one of the new html composite data types for Rgenetics
        Dan Blankenberg with mods by Ross Lazarus 
        October 2007
    """
    out = file(outfile,'w')
    out.write(galhtmlprefix % progname)

    if len(flist) > 0:
        out.write('<ol>\n')
        for i, data in enumerate( flist ):
           out.write('<li><a href="%s" mime-type="application/binary">%s</a></li>\n' % \
                     (os.path.split(data)[-1],os.path.split(data)[-1]))
        out.write('</ol>\n')
    else:
           out.write('No files found')
    out.write('<hr>Log:<br>')
    out.write('\n<br>'.join(log))
    out.write('<hr>\n')
    out.write("</div></body></html>")
    out.close()


       
def snpMatrixcon():
    """
         <command interpreter="python">
            rgSNPMatrixMaker.py '$input_file.extra_files_path' '$input_file.metadata.base_name'
            '$outfile1' '$outfile1.files_path' '$title' '$logpath'
        </command>
        
    """
    inpath, basename,outfile,outfilepath,title,logpath,rexe = sys.argv[1:7]
    infile = os.path.join(inpath, basename)
    tempMapName = '%s.info' % infile
    pedfile = '%s.ped' % infile
    mapname = '%s.map' % infile
    smname = '%s.snpmatrix' % title
    log = []
    try:
        os.makedirs(outfilepath)
    except:
        pass
    outsm = os.path.join(outfilepath,smname)
    if os.path.exists(pedfile):
        try:
            os.makedirs(outfilepath)
        except:
            pass
        ss = '%s%s' % (string.punctuation,string.whitespace)
        ptran =  string.maketrans(ss,'_'*len(ss))
        title = title.translate(ptran)
        log.append('# %s got %s\n' % (progname,sys.argv))
        if not os.path.exists(tempMapName): # make it
            if os.path.exists(mapname):
                log.append('Making info file for snpmatrix %s\n' % (tempMapName))
                mf = file(mapname,'r')
                m = mf.readlines()
                markers = [x.strip().split() for x in m] # chrom snp gen abs
                tempMap = file(tempMapName,'w')
                map = ['%s %s %s\n' % (x[1],x[3],x[0]) for x in markers]
                # snp,abspos,chrom in genomic order for haploview
                tempMap.write(''.join(map))
                tempMap.close()
            else:
                log.append('## cannot find %s - so cannot make %s - Sorry\n' % (mapname,tempMapName))
                sys.exit(1)
        shandle,sname = tempfile.mkstemp('.R','script')     
        ehandle,ename = tempfile.mkstemp('.log','elog') # capture stderr!
        ohandle,oname = tempfile.mkstemp('.log','log')
        tempo = file(oname,'w')
        tempe = file(ename,'w')
        f = file(sname,'w')
        s = 'library("snpMatrix")\n'
        f.write(s)
        s = '%s = read.snps.pedfile(file="%s")\n' % (title,pedfile)
        f.write(s)
        s = 'save(%s, file="%s")\n' % (title,outsm)
        f.write(s)
        f.close()
        vcl = [rexe,' --vanilla <',sname] #,' >',oname]
        x = subprocess.Popen(' '.join(vcl),stdout=tempo,stderr=tempe,shell=True)
        retval = x.wait()
        tempo.close()
        tempe.close()
        log.append('stdout from R:\n')
        log += file(oname,'r').readlines()
        log.append('stderr from R:\n')
        log += file(ename,'r').readlines()
        log.append('\n##retval = %s\n' % (str(retval)))
        os.unlink(sname)
        os.unlink(oname)
        doImport(outfile,[smname],log)
        log.append('Imported %s into %s\n' % (smname,outfile))
        logf = file(logpath,'w')
        logf.write('\n'.join(log))
        logf.close()
    else:
        print >> sys.stdout, '## Cannot open inputfilepath %s - exists and readable by you?' % pedfile
        print >> sys.stdout, '## Usage: %s inpath basename outfile outfilepath title logpath' % progname
 

if __name__ == "__main__":
    print >> sys.stdout, '# rgSNPMatrixMaker got %s\n' % sys.argv
    if len(sys.argv) >= 7:
        snpMatrixcon()           
    else:
        print >> sys.stdout, '## Missing one or more command line parameters'
        print >> sys.stdout, '## Usage: %s inpath basename outfile outfilepath title logpath' % progname

