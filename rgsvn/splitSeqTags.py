#!/usr/bin/python
## hacked from Michael Cho's script 
## mipBarcodes.py
## examines each line, then creates individual files for each barcode.

import string, sys, os, gzip

pcpgmbcl = ['ACAGGC', 'ACATCA', 'ACGAGG', 'ACTACT', 'AGCTGA', 'AGGCAG', 'AGGTCT', 'AGTGAC', 'ATATTG', 
'ATCACG', 'ATGGAT', 'ATTCGA', 'CACTCA', 'CAGAGT', 'CAGCTG', 'CATGGA', 'CGGCGC', 'CGGTTA', 'CGTCAT', 
'CGTGTG', 'CTCATT', 'CTGGCA', 'CTGTAG', 'CTTAGG', 'GACAGG', 'GAGGAC', 'GAGTCG', 'GATTAT', 'GCGGCT',
 'GCGTTC', 'GCTGTA', 'GCTTGG', 'GTGGTG', 'GTGTGT', 'GTTCTT', 'GTTGGC', 'NNNNNN', 'TACACT', 'TACTGC',
 'TAGGTA', 'TATCAG', 'TCACTG', 'TCATGT', 'TCGGAG', 'TCTCGC', 'TGGCTT', 'TGGTGG', 'TGTGGT', 'TGTTCA']


bcf = {}

def splitSTags(finname=None,foutname=None, bcl=pcpgmbcl):
    """
    """
    assert os.path.isfile(finname)
    bcf = [open('%s_%s.fastq' % (fileRoot,x),'w') for x in bcl] # create a lookup dict of files
    bc = dict(zip(bcl,bcf)) 
    bcCounts = {}
    print os.path.splitext(finname)
    if os.path.splitext(finname)[-1] in ['.gz',]:
        fin = gzip.open(finname, 'r')
    else:
        fin = open(finname, 'r')
    for (i, line) in enumerate(fin):
        if i % 4 == 0:
        ## if (line[0]=="@"):
            barcode=line.split("#")[1][:6]
            ## first count the barcodes
            bcCounts.setdefault(barcode,0)
            bcCounts[barcode]+=1
            ## assign it a dummy value if it is not in official list
            if not(bc.has_key(barcode)):
                barcode='NNNNNN'
        if (i+1) % 1000000 == 0:
            print 'On line',(i+1)
        bc[barcode].write(line)
    for f in bc.itervalues():
        f.close()

    ## file("barcodeCounts.txt","w").write(repr(bcCounts))
    outlog = open(fileRoot+"BarcodeCounts.xls","w")
    outlog.write('SeqTag\tCount\n')
    for key in sorted(bcCounts):
        n = bcCounts[key]
        outlog.write('%s\t%d\n' % (str(key),n))
        if bc.get(key,None): # is a target
           print '%s\%d' % (key,bcCounts[key])
    outlog.close()


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print '## Need an input file to be split by sequence tag'
        sys.exit(1)
    fin=sys.argv[1]

    if len(sys.argv) > 2:
        fileRoot=sys.argv[2]
    else:
        fileRoot=sys.argv[1]
        fileRoot=fileRoot.split(".")[0]
    fout = "%s_BarcodeCounts.txt" % fileRoot
    splitSTags(finname=fin,foutname=fout,bcl=pcpgmbcl)

