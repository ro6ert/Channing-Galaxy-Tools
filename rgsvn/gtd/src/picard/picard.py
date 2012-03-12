# utilities for picard
#
# copyright 2010 ross lazarus
# released under the LGPL
#

import os
import sys
import datetime
import string
from subprocess import Popen

from mako import Template


class Picard(object):
    JAR_DIR = os.path.join(
        os.getenv('GALAXY_DATA_INDEX_DIR'), 'shared', 'jars', 'CalculateHsMetrics.jar')

    def __init__(self, input_file, **kw):
        self.input = input_file
        self.output_file = kw.pop('output_file', None)
        self.output_dir = kw.pop('output_dir', None)
        self.bait = kw.pop('bait_file', None)
        self.target = kw.pop('target_file', None)



    def hs_metrics(self, job_name):
        killme = string.punctuation + string.whitespace
        trantab = string.maketrans(killme,'_'*len(killme))
        title = job_name.translate(trantab)

        
    def make_HS_html_report(self):
        pass

    def collect_GC_bias(self):
        pass

    def AS_stats(self):
        pass



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
galhtmlattr = """<h3><a href="http://rgenetics.org">Rgenetics</a> tool %s run at %s</h3>"""
galhtmlpostfix = """</div></body></html>\n"""



def fixPicardOutputs(tempout=None,output_dir=None,log_file=None,html_output=None,progname=None,cl=[],transpose=True):
    """
    picard produces long hard to read tab header files
    make them available but present them transposed for readability
    """
    rstyle="""<style type="text/css">
    tr.d0 td {background-color: oldlace; color: black;}
    tr.d1 td {background-color: aliceblue; color: black;}
    </style>"""    
    cruft = []
    dat = []    
    try:
        r = open(tempout,'r').readlines()
    except:
        r = []
    for row in r:
        if row.strip() > '':
            srow = row.split('\t')
            if row[0] == '#':
                cruft.append(row.strip()) # want strings
            else:
                dat.append(srow) # want lists
    
    res = [rstyle,]
    res.append(galhtmlprefix % progname)   
    res.append(galhtmlattr % (progname, datetime.datetime.now()))
    flist = os.listdir(output_dir) 
    res.append('<b>Your job produced the following output files.</b><hr/>\n')
    res.append('<table>\n')
    for i,f in enumerate(flist):
         fn = os.path.split(f)[-1]
         res.append('<tr><td><a href="%s">%s</a></td></tr>\n' % (fn,fn))
    res.append('</table><p/>\n') 
    res.append('<b>Picard log</b><hr/>\n') 
    try:
        rlog = ['<pre>',]
        rlog += open(log_file,'r').readlines()
        rlog.append('</pre>')
        res += rlog
    except:
        res.append("Odd, Picard left no log file %s - must have really barfed badly?" % log_file)
    if len(cruft) + len(dat) > 0:
        if transpose:
            res.append('<b>Picard output (transposed for readability)</b><hr/>\n')       
        else:
            res.append('<b>Picard output</b><hr/>\n')       
        res.append('<table cellpadding="3" >\n')
        if len(cruft) > 0:
            cres = ['<tr class="d%d"><td>%s</td></tr>' % (i % 2,x) for i,x in enumerate(cruft)]
            res += cres
        if len(dat) > 0: 
            maxrows = 100
            if transpose:
                tdat = map(None,*dat) # transpose an arbitrary list of lists
                missing = len(tdat) - maxrows
                tdat = ['<tr class="d%d"><td>%s</td><td>%s</td></tr>\n' % ((i+len(cruft)) % 2,x[0],x[1]) for i,x in enumerate(tdat) if i < maxrows] 
                if len(tdat) > maxrows:
                   tdat.append('<tr><td colspan="2">...WARNING: %d rows deleted for sanity...see raw files for all rows</td></tr>' % missing)
            else:
                tdat = ['<tr class="d%d"><td>%s</td></tr>\n' % ((i+len(cruft)) % 2,x) for i,x in enumerate(dat) if i < maxrows] 
                if len(dat) > maxrows:
                    missing = len(dat) - maxrows      
                    tdat.append('<tr><td>...WARNING: %d rows deleted for sanity...see raw files for all rows</td></tr>' % missing)
            res += tdat
        res.append('</table>\n')   
    else:
        res.append('<b>No Picard output found - please consult the Picard log above for an explanation</b>')
    res.append('<hr/>The freely available <a href="http://picard.sourceforge.net/command-line-overview.shtml">Picard software</a> \n') 
    res.append( 'generated all outputs reported here, using this command line:<br/>\n<pre>%s</pre>\n' % ''.join(cl))   
    res.append(galhtmlpostfix) 
    outf = open(html_output,'w')
    outf.write(''.join(res))   
    outf.write('\n')
    outf.close()

def bedToPicInterval(infile=None,outfile=None):
    """
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
    """
    # bedToPicard.py
    # ross lazarus October 2010
    # LGPL 
    # for Rgenetics
    """
    collect header info and rewrite bed with header for picard
    """
    featlen,fk,otherHeaders = getFlen(bedfname=infile)
    try:
        outf = open(outfile,'w')
    except:
        print '###ERROR: writePic unable to open output picard file %s' % outfile
        sys.exit(1)
    inf = open(infile,'r') # already tested in getFlen
    header = ['@SQ\tSN:%s\tLN:%d' % (x,featlen[x]) for x in fk]
    if len(otherHeaders) > 0:
        header += otherHeaders
    outf.write('\n'.join(header))
    outf.write('\n')
    for row in inf:
        row = row.strip()
        if len(row) > 0: # convert zero based start coordinate to 1 based
            if row.startswith('@'):
                continue
            else:
                srow = row.split('\t')
                srow[1] = '%d' % (int(srow[1])+1) # see http://genome.ucsc.edu/FAQ/FAQtracks.html#tracks1
                outf.write('\t'.join(srow))
                outf.write('\n')
    outf.close()
    inf.close()

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

def getFlen(bedfname=None):
    """
    find all features in a BED file and sum their lengths
    """
    features = {}
    otherHeaders = []
    try:
        infile = open(bedfname,'r')
    except:
        raise IOError, 'getFlen unable to open bedfile %s' % bedfname
    for i,row in enumerate(infile):
        if row.startswith('@'): # add to headers if not @SQ
            if not row.startswith('@SQ'):
                otherHeaders.append(row)
        else:
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
    fk = features.keys()
    fk = sorted(fk, key=keynat)
    return features,fk,otherHeaders
