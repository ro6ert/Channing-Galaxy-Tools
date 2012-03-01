#!/usr/bin/env python
# This code exists in 2 places: ~/datatypes/converters and ~/tools/filters
import sys,os,tempfile,subprocess,glob,shutil

assert sys.version_info[:2] >= ( 2, 4 )

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

def doPbedToLpedConvert(destdir=None,sourcepath=None,basename=None,logf=None):
    """
    """
    outroot = os.path.join(destdir,basename)
    missval = getMissval(inped = pedf)
    cl = 'plink --noweb --bfile %s --recode --out %s > %s' % (sourcepath,outroot,missval,logf)
    p = subprocess.Popen(cl,shell=True,cwd=destdir)
    retval = p.wait() # run plink
    return 
   
def doConvert(destdir=None,sourcepath=None,basename=None,logf=None):
    """
    """
    pedf = '%s.ped' % sourcepath
    outroot = os.path.join(destdir,basename)
    missval = getMissval(inped = pedf)
    if not missval:
        print '### lped_to_pbed_converter.py cannot identify missing value in %s' % pedf
        missval = '0'
    cl = 'plink --noweb --file %s --make-bed --out %s --missing-genotype %s > %s' % (sourcepath,outroot,missval,logf)
    p = subprocess.Popen(cl,shell=True,cwd=destdir)
    retval = p.wait() # run plink
    return missval

def __main__():
    progname = os.path.split(sys.argv[0])[-1]
    input_name = sys.argv[1]
    output_name = sys.argv[2]
    output_path = sys.argv[3]
    try:
        os.makedirs(output_path)
    except:
        pass # no need
    input_datatype = sys.argv[4] # eg lped, pbed
    output_datatypes = sys.argv[5].split(',') # 
    basefile = os.path.split(input_name)[-1]
    basename = os.path.splitext(basefile)[0]
    tfd,tlog = tempfile.mkstemp(prefix='rgLpedPbedlog')
    if input_datatype='lped': # linkage ped in
        mappath = '%s.map' % os.path.splitext(input_name)[0]
        mapcopy = os.path.join(output_path,'%s.map' % basename)
        # need this for caco
        shutil.copy(mappath,mapcopy)
    if 'pbed' in output_datatypes:
        missval = doConvert(destdir=output_path,sourcepath=input_name,basename=basename,logf=tlog)
        plinklog = file(tlog,'r').read()
        os.unlink(tlog)
        info = 'Converted %s, using missing value code = %s\n' % (basename,missval)
    if 'lped' in output_datatypes:
        doPbedToLpedConvert(destdir=output_path,sourcepath=input_name,basename=basename,logf=tlog)
        plinklog = file(tlog,'r').read()
        os.unlink(tlog)
        info = 'Converted %s, using missing value code = %s\n' % (basename,missval)
        
    out = file(output_name,'w') # start again
    out.write(galhtmlprefix % progname)
    out.write('<h2>Rgenetics data imported from %s</h2>\n<hr><ul>' % basename)
    flist = glob.glob(os.path.join(output_path, '%s.*' % (basename))) # all the files matching basename
    for i, data in enumerate( flist ):
        ofn = os.path.split(data)[-1]
        out.write('<li><a href="%s" type="application/binary">%s</a></li>\n' % (ofn,ofn))
    out.write('<hr><h3>Plink log follows:</h3><pre>')
    out.write(plinklog)
    out.write('</pre><br>')
    out.write("</div></body></html>")
    out.close()
    print info 

if __name__ == "__main__": __main__()
