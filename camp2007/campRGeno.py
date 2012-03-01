"""
released under the terms of the LGPL
copyright ross lazarus August 2007 
for the rgenetics project

Special galaxy tool for the camp2007 data
Allows grabbing genotypes from an arbitrary region

Needs a mongo results file in the location hardwired below or could be passed in as
a library parameter - but this file must have a very specific structure
rs chrom offset float1...floatn

called as

  
    <command interpreter="python">
        campRGeno.py $region "$rslist" "$title" $output1 $log_file $userId
    </command>

"""


import sys, array, os

progname = os.path.split(sys.argv[0])[1]

DATA_FILE = '/usr/local/galaxy/data/camp550KV3/550KV3_PBAT_auto.ped'
MAP_FILE = '/usr/local/galaxy/data/camp550KV3/camp2007.map'
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


def doImport(outfile='test',flist=[]):
    """ import into one of the new html composite data types for Rgenetics
        Dan Blankenberg with mods by Ross Lazarus 
        October 2007
    """
    out = open(outfile,'w')
    out.write(galhtmlprefix % progname)

    if len(flist) > 0:
        out.write('<ol>\n')
        for i, data in enumerate( flist ):
           out.write('<li><a href="%s">%s</a></li>\n' % (os.path.split(data)[-1],os.path.split(data)[-1]))
        out.write('</ol>\n')
    else:
           out.write('No files found')
    out.write("</div></body></html>")
    out.close()


def subset():
    """  ### Sanity check the arguments
        <command interpreter="python">
            campRGeno.py $region "$rslist" "$title" $output1 $log_file $userId
        </command>
    """

    if len(sys.argv) < 6: 
        print >> sys.stdout, '##!campRGeno: Expected 7 params in sys.argv, got %d (%s)' % (len(sys.argv), sys.argv)
        sys.exit(1)

    ### Figure out what genomic region we are interested in
    region = sys.argv[1]
    orslist = sys.argv[2].replace('X',' ').lower() # galaxy replaces newlines with XX - go figure
    title = sys.argv[3].replace(' ','_') # for outputs
    outfile = sys.argv[4]
    log_file = sys.argv[5]
    userId = sys.argv[6] # for library
    lf = file(log_file,'w')

    if region > '':
        try: # TODO make a regexp?
            c,rest = region.split(':')
            c = c.replace('chr','')
            rest = rest.replace(',','') # remove commas
            spos,epos = rest.split('-')
            spos = int(spos)
            epos = int(epos)
            print >> sys.stdout, '## %s parsing chrom %s from %d to %d' % (progname,c,spos,epos)
        except:
            print >> sys.stdout, '##! %s unable to parse region %s - MUST look like "chr8:10,000-100,000' % (progname,region)
            sys.exit(1)
        ### Figure out which markers are in this region
        mfile = open(MAP_FILE, 'r')
        markers = []
        for line in mfile:
            line = line.strip()
            if not line: continue
            chrom, snp, genpos, abspos = line.split()
            abspos = int(abspos)
            if chrom == c and (spos <= abspos <= epos):
                markers.append((chrom,abspos,snp)) # decorate for sort into genomic
        markers.sort()
        rslist = [x[2] for x in markers] # drop decoration
        rsdict = dict(zip(rslist,rslist))
        mfile.close()
    else:
        rslist = orslist.split() # galaxy replaces newlines with XX - go figure
    if len(rslist) == 0:
            print >> sys.stdout, '##! %s found no rs numbers in' % (progname,sys.argv[1:3])
            sys.exit(1)        
    rsdict = dict(zip(rslist,rslist)) # for quick lookups
    s = '## %s looking for %d rs (%s)' % (progname,len(rslist),rslist[:5])
    print >> sys.stdout, s
    ### Now parse the first line of the ped file to determine which columns
    ### have the markers we want ...
    try:
        dfile = open(DATA_FILE, 'r')
    except: # bad input file name?
        print >> sys.stdout, '##!campRG unable to open file %s' % (DATA_FILE)
        raise
        sys.exit(1)
    header = dfile.readline().strip()
    snps = header.split()
    snpcols = {}
    for column, snp in enumerate(snps):
        if rsdict.get(snp,None):
            snpcols[snp] = column
    wewant = [(6+(2*snpcols[x])) for x in rslist] # 
    # column indices of first geno of each marker pair to get the markers into genomic
    ### ... and then parse the rest of the ped file to pull out
    ### the genotypes for all subjects for those markers
    # /usr/local/galaxy/data/rg/1/lped/
    outext = 'lped'
    libroot = os.path.join(library,userId)
    outpath = os.path.join(libroot,outext)
    try:
        os.makedirs(outpath)
    except:
        pass
    outname = os.path.join(outpath,title)
    mapfname = '%s.map' % outname
    pedfname = '%s.ped' % outname
    ofile = file(pedfname, 'w')

    # make a map file in the lped library
    mf = file(mapfname,'w')
    map = ['%s\t%s\t0\t%d' % (x[0],x[2],x[1]) for x in markers] # chrom,abspos,snp in genomic order
    mf.write('%s\n' % '\n'.join(map))
    mf.close()
    nrows = 0
    for line in dfile:
        line = line.strip()
        if not line:
            continue
        fields = line.split()
        preamble = fields[:6]
        g = ['%s %s' % (fields[snpcol], fields[snpcol+1]) for snpcol in wewant]
        ofile.write('%s %s\n' % (' '.join(preamble), ' '.join(g)))
        nrows += 1
    dfile.close()    
    ofile.close()
    s = '## %s: wrote %d markers, %d subjects for region %s\n' % (progname,len(rslist),nrows,region)
    lf.write(s)
    print >> sys.stdout,s
    lf.close()
    doImport(outfile,[mapfname,pedfname])

if __name__ == "__main__":
    subset()

