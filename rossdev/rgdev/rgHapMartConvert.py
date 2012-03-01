"""
changed to create separate haploview 'hapmap' format 
files for each chromosome so haploview doesn't barf

Note - currently only works for haploview beta release
as hapmap 3 data output format was used as the model...

Ross Lazarus 30 June 2009
Hapmart dump to haploview hapmap genotype format for Benji

logic is:

Just to be clear - check this over please?

Hapmap format

1:rs# 2:alleles 3:chrom 4:pos 5:strand 6:assembly# 7:center 8:protLSID
9:assayLSID 10:panelLSID 11:QCcode 12:NA06984 NA0698..

HapMart format
1:chromosome   2:position    3:Strand   4:marker id    5:Alleles
6:reference allele    7:genotyping center    8:genotyping platform
9:POP:CEU 10[NA06984...

So, from hapmart to hapmap
1->3
2->4
3->5
4->1
5->2
     6 is junk
     7 is junk
     8 is junk
     9 is junk
     10 is junk
     11 is QC+
10->12 with tabs as seps instead of space and removing the [] chars
"""

import sys,time,os

PROGNAME = os.path.split(sys.argv[0])[-1]

def timenow():
    """return current time as a string
    """
    return time.strftime('%d/%m/%Y %H:%M:%S', time.localtime(time.time()))

def writeHapmap(outfpath=None,res=[], infpath=None,hmhead=None):
    """
    rs# alleles chrom pos strand assembly# center protLSID assayLSID panelLSID QCcode NA06984 NA06985 NA0$

    hapmap format looks like:
    #Tue Jun 30 11:15:24 2009: HapMap genotype data dump, SNPs genotyped in population CEU on chr4:38450658..38460984
#For details on file format, see http://www.hapmap.org/genotypes/
rs# alleles chrom pos strand assembly# center protLSID assayLSID panelLSID QCcode NA06984 NA06985 NA06986 NA06989 NA06991 NA06993 NA06994 NA06995 NA06997 NA07000 NA07014 NA07019 NA07022 NA07029 NA07031 NA07034 NA07037 NA07045 NA07048 NA07051 NA07055 NA07056 NA07345 NA07346 NA07347 NA07348 NA07349 NA07357 NA07435 NA10830 NA10831 NA10835 NA10836 NA10837 NA10838 NA10839 NA10840 NA10843 NA10845 NA10846 NA10847 NA10850 NA10851 NA10852 NA10853 NA10854 NA10855 NA10856 NA10857 NA10859 NA10860 NA10861 NA10863 NA10864 NA10865 NA11829 NA11830 NA11831 NA11832 NA11839 NA11840 NA11843 NA11881 NA11882 NA11891 NA11892 NA11893 NA11894 NA11917 NA11918 NA11919 NA11920 NA11930 NA11931 NA11992 NA11993 NA11994 NA11995 NA12003 NA12004 NA12005 NA12006 NA12043 NA12044 NA12045 NA12056 NA12057 NA12144 NA12145 NA12146 NA12154 NA12155 NA12156 NA12234 NA12236 NA12239 NA12248 NA12249 NA12264 NA12272 NA12273 NA12275 NA12282 NA12283 NA12286 NA12287 NA12335 NA12336 NA12340 NA12341 NA12342 NA12343 NA12344 NA12347 NA12348 NA12375 NA12376 NA12383 NA12386 NA12399 NA12400 NA12413 NA12489 NA12546 NA12707 NA12708 NA12716 NA12717 NA12718 NA12739 NA12740 NA12748 NA12749 NA12750 NA12751 NA12752 NA12753 NA12760 NA12761 NA12762 NA12763 NA12766 NA12767 NA12775 NA12776 NA12777 NA12778 NA12801 NA12802 NA12812 NA12813 NA12814 NA12815 NA12817 NA12818 NA12827 NA12828 NA12829 NA12830 NA12832 NA12842 NA12843 NA12864 NA12865 NA12872 NA12873 NA12874 NA12875 NA12877 NA12878 NA12889 NA12890 NA12891 NA12892
rs10776482 A/G chr4 38451180 + ncbi_b36 sanger urn:LSID:illumina.hapmap.org:Protocol:Human_1M_BeadChip:3 urn:LSID:sanger.hapmap.org:Assay:H1Mrs10776482:3 urn:lsid:dcc.hapmap.org:Panel:CEPH-60-trios:3 QC+ AA AG AA AA AG AA AA AG AA AA AA NN AG AA AA NN AG AA NN AA AA AG AG AG AA AG AA AG AA AG AA AG AA AA AA AA AG AA AA AG AA AA NN AA AA AG AA AA NN AA NN AA AA AG AA AG AG AA AG AG AA AA AA AG AA AA AA AG AG AA AA AA AA AA AA AA AA AA AA NN AA AA AA AG AA AA AA AG AG AA AA AA AG AG NN AG GG AG AA AA AA AA AG AG AG AA AA AG AA AA AG GG AA AA AA AG AA GG AA AA AG AG AA AG AG AA AG NN AA AA AG AG AA AA GG AG AG GG AA AG AG AA AA AA AA AG AG AG AG AG AA AA GG AA AA AA AA AG AA AA AA AA AA AG AA AG AG AG AA AG AA AG GG AA
rs4129009 C/T chr4 38451284 + ncbi_b36 sanger urn:LSID:illumina.hapmap.org:Protocol:Human_1M_BeadChip:3 urn:LSID:sanger.hapmap.org:Assay:H1Mrs4129009:3 urn:lsid:dcc.hapmap.org:Panel:CEPH-60-trios:3 QC+ TT CT TT TT CT TT 
    """
    head = ' '.join(hmhead)
    f = open(outfpath,'w')
    t = ['# %s generated by %s from infile %s' % (timenow(),PROGNAME,infpath)]
    t.append('# see whatever')
    f.write('\n'.join(t))
    f.write('\n')
    f.write(head)
    f.write('\n')
    f.write('\n'.join(res))
    f.write('\n')
    f.close()

def readMart(infpath=None):
    """ read header row and parse what we need
    hapmart format looks like:
    chromosome    position    Strand    marker id    Alleles    reference allele    genotyping center    genotyping platform    POP:CEU [NA06984 NA06985 NA06986 NA06989 NA06991 NA06993 NA06994 NA06995 NA06997 NA07000 NA07014 NA07019 NA07022 NA07029 NA07031 NA07034 NA07037 NA07045 NA07048 NA07051 NA07055 NA07056 NA07345 NA07346 NA07347 NA07348 NA07349 NA07357 NA07435 NA10830 NA10831 NA10835 NA10836 NA10837 NA10838 NA10839 NA10840 NA10843 NA10845 NA10846 NA10847 NA10850 NA10851 NA10852 NA10853 NA10854 NA10855 NA10856 NA10857 NA10859 NA10860 NA10861 NA10863 NA10864 NA10865 NA11829 NA11830 NA11831 NA11832 NA11839 NA11840 NA11843 NA11881 NA11882 NA11891 NA11892 NA11893 NA11894 NA11917 NA11918 NA11919 NA11920 NA11930 NA11931 NA11992 NA11993 NA11994 NA11995 NA12003 NA12004 NA12005 NA12006 NA12043 NA12044 NA12045 NA12056 NA12057 NA12144 NA12145 NA12146 NA12154 NA12155 NA12156 NA12234 NA12236 NA12239 NA12248 NA12249 NA12264 NA12272 NA12273 NA12275 NA12282 NA12283 NA12286 NA12287 NA12335 NA12336 NA12340 NA12341 NA12342 NA12343 NA12344 NA12347 NA12348 NA12375 NA12376 NA12383 NA12386 NA12399 NA12400 NA12413 NA12489 NA12546 NA12707 NA12708 NA12716 NA12717 NA12718 NA12739 NA12740 NA12748 NA12749 NA12750 NA12751 NA12752 NA12753 NA12760 NA12761 NA12762 NA12763 NA12766 NA12767 NA12775 NA12776 NA12777 NA12778 NA12801 NA12802 NA12812 NA12813 NA12814 NA12815 NA12817 NA12818 NA12827 NA12828 NA12829 NA12830 NA12832 NA12842 NA12843 NA12864 NA12865 NA12872 NA12873 NA12874 NA12875 NA12877 NA12878 NA12889 NA12890 NA12891 NA12892]
chr2    230747154    +    rs1004868    A/G    G    sanger    Illumina_1M    GG NN AG GG GG GG GG GG GG AG GG NN GG GG GG NN GG GG NN GG GG GG GG GG AG GG AG GG GG GG GG GG GG GG AG GG GG AG GG AG GG GG NN GG GG GG GG GG NN GG NN GG GG GG GG GG GG GG GG GG GG GG GG GG GG GG GG GG GG GG AG GG GG GG GG GG AG GG GG NN GG GG GG GG GG GG GG GG AG GG AG GG GG GG NN GG GG GG AG GG AG AG GG GG GG GG GG AG GG GG AG AG GG GG GG AG GG GG GG GG GG GG GG GG GG GG AG NN AG GG GG GG GG GG AG GG GG GG GG GG GG GG AG GG GG AG GG GG GG AG GG GG GG GG GG GG GG GG GG AG GG AG GG GG GG GG GG GG GG AG AG GG AG AG
chr2    230747515
move 1->3
2->4
3->5
4->1
5->2
    """
    print 'infpath=',infpath
    f = open(infpath,'r')    
    filler = ['ncbi_b36','sanger','urn:LSID:illumina.hapmap.org:Protocol:Human_1M_BeadChip:3',
              'urn:LSID:sanger.hapmap.org:Assay:H1Mrs10776482:3','urn:lsid:dcc.hapmap.org:Panel:CEPH-60-trios:3','QC+']
    res = {}
    for i,row in enumerate(f):
        if i == 0:
            header = row.strip().split('\t') # spaces in ids dealt with later
            idlist = header[-1].strip()
            h = [x.lower() for x in header] # mmmff
            chrpos = h.index('chromosome')
            pospos = h.index('position')
            strandpos = h.index('strand')
            rspos = h.index('marker id')
            allelepos = h.index('alleles')
            wewant = [rspos,allelepos,chrpos,pospos,strandpos]
            idlist = idlist.replace('[','')
            idlist = idlist.replace(']','')
            idlist = idlist.split()[1:] #get rid of pop:Ceu
            nid = len(idlist)
            hmhead = 'rs# alleles chrom pos strand assembly# center protLSID assayLSID panelLSID QCcode'.split()
            hmhead += idlist           
        else:
            rowl = row.strip().split()
            if len(rowl) > 2:
                genos = rowl[-nid:] # last nid things
                hrow = [rowl[x] for x in wewant] # gets the fields in order
                chrom = rowl[chrpos]
                hrow += filler
                hrow += genos
                res.setdefault(chrom,[])
                res[chrom].append(' '.join(hrow))
    return res,hmhead
    
    
if __name__ == "__main__":   
    if len(sys.argv) < 2:
        print 'Need a hapmart output file as the first command line parameter'
        sys.exit(1)
    infname = sys.argv[1]
    fname = os.path.split(infname)[-1]
    outfname = 'hapmap_%s' % fname
    res,hmhead = readMart(infpath=infname)
    chroms = res.keys()
    chroms.sort()
    for c in chroms:
        outf = '%s_%s' % (c,fname)
        writeHapmap(outfpath=outf,res=res[c],infpath=infname,hmhead=hmhead)
    
    
    