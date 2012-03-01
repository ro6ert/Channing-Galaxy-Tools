"""
TODO: Note we need to add proper blat annotation to the affy annotation already here.

We currently convert xY (eg TY) to xx (eg TT). 

This probably only sometimes makes sense, so you should think about how you want to handle 
X and X-psuedo autosomal in particular. Any suggestions welcomed.

NEW version for GAIN control munging
ross lazarus feb 2 2009
note files are all gz still - we trade off cpu cycles for disk i/o using
a gz file reader. TODO - test to see how this trade off works?

NO pedigree file available - instead, use the GAINcontrol_list.csv supplied by Brooke

SubjID	individual_id	site_number	sex	sampling_age	black	white	control												
738	11365	150	1	65	1	0	1															
793	11467	150	1	68	0	1	1															
853	11565	150	1	48	0	1	1															

------------------------------------------------------------------
dbgap (framingham share, Gain controls are the two I've seen so far)
data arrives as a mess of gzipped individual genotype files
inside tarballs inside RAR archives.
Oy

Matrix format has less redundancy than individual format
but you need individual for intensities

Code here reads indiv format - only - old matrix code cleaned out for clarity. 

The individual data contains
intensity and quality scores that are useful for evaluating hits, but the scale makes managing
these challenging. eg sorting the intensity data for only 1/3 (the nonprofit release) takes about 10
hours. Not sure whether it's worth the bother of loading it up...

"""

import glob, sys, os, MySQLdb, csv, copy, array, gzip, time
debug = 0 # prints mysterious rs numbers

def readControls(fname='GAINcontrol_list.csv'):
    """Brooke provided this csv file in an email, late January 2008
       we extract race etc for the samples we want - use the list to control
       the ped file builder


==> GAINcontrol_list.csv <==
SubjID,individual_id,site_number,sex,sampling_age,black,white,control
738,11365,150,1,65,1,0,1
793,11467,150,1,68,0,1,1
853,11565,150,1,48,0,1,1
857,11571,150,1,46,1,0,1
878,11670,150,1,40,1,0,1
925,11762,150,1,71,0,1,1
997,11878,150,2,34,0,1,1
1033,11941,150,2,62,1,0,1
1048,11966,150,1,45,0,1,1


    """
    controls = {} # keep id's for subjects keyed here
    cped = {} # keep fake ped here
    dkeys = ['SubjID','iid','site','sex','age','black','white','control']
    f = file(fname,'r')
    fl = f.readlines()
    head = fl.pop(0).strip().split(',')
    for nl in fl:
        ll = nl.strip().split(',')
        if len(ll) == 8:
            SubjID,individual_id,site_number,sex,sampling_age,black,white,control = ll
            if control == '1':
                d = dict(zip(dkeys,ll))
                controls[individual_id] = d
                fakeped = [individual_id,individual_id,'0','0',sex,'1'] # all unaffected - hope the genders are right
                cped[individual_id] = fakeped
        else:
          print '## line %d len %d - looking for 8' % (lnum,len(ll))   
    print 'readControls found',len(cped)
    return controls,cped


def makeAffyMap(annoroot='./', gdir='./', outfpref=None):
    """
    eeesh. duplicate fucking markers in the output.
    hacked 24 march rml
    john repaired this
    refactored so drop in replacement for getMarkers

    Leaving this here in case it's easy to fix
    I left this cos I couldn't get it to work quickly - so I parsed the anno file

    No idea how the affy id's are arranged and quicker to just
    use the crap annotation file directly
    need a map file - john has blatted the affy 6 probes and created 2 files copied from meme
'
I knew I had done this before:

meme:/home/SHARP/IMPORTER/GenomeWideSNP_6_flanking_sequences.map

and in that same directory is a file which maps
from affy's special names to rs numbers:

affyToRSMap.txt
'

[rerla@beast dbgap]$ head affyToRSMap.txt 
AFFX-SNP_9998494 rs1075681
AFFX-SNP_9997271 rs7913189
AFFX-SNP_9997260 rs6831487
AFFX-SNP_9986312 rs869487
AFFX-SNP_9985333 rs6519955
AFFX-SNP_998040 rs444163
AFFX-SNP_9976462 rs2810434
AFFX-SNP_9976431 rs9304712
AFFX-SNP_9963677 rs6029111
AFFX-SNP_99759 rs479341
[rerla@beast dbgap]$ head GenomeWideSNP_6_flanking_sequences.map 
-1	SNP_A-1808348	0	-1
-1	SNP_A-1895632	0	-1
-1	SNP_A-1904584	0	-1
-1	SNP_A-1908397	0	-1
-1	SNP_A-1910995	0	-1
-1	SNP_A-2155895	0	-1
-1	SNP_A-2156784	0	-1
-1	SNP_A-4241806	0	-1
-1	SNP_A-8291010	0	-1
-1	SNP_A-8425903	0	-1

    """
    f = file('affyToRSMap.txt','r').readlines()
    fl = [x.strip().split() for x in f]
    fl = [x for x in fl if len(x) >= 2]
    #a = ['SNP_A-%s' % x[0].split('_')[-1] for x in fl] # get rid of bogus AFFX-SNP to replace with bogus bogus SNP_A-
    a = [x[0] for x in fl] # get rid of bogus AFFX-SNP to replace with bogus bogus SNP_A-
    #print 'a=',a[:20]
    rs = [x[1] for x in fl]
    ators = dict(zip(a,rs))

    f = file('GenomeWideSNP_6_flanking_sequences.map','r').readlines()
    fl = [x.split() for x in f]
    fl = [x for x in fl if len(x) > 0 and int(x[0])>0] # ignore missing chroms?
    print 'fl=',fl[:20]
    affyToMap = {}
    rslist = [] # for decorated chrom offset sorting
    rsseen = {} # for stoopid dupe rs
    notfound = 0
    for x in fl:
        chrom,aid,gpos,ppos = x
        rs = ators.get(aid, None)
        if rs:
             if rsseen.get(rs,None):
                # damn duplicates in John's inputs - ? multiple probes 
                print '# ignoring second %s, rs=%s' % (aid,rs)
             else:
                rsseen[rs] = rs
                affyToMap[aid] = [chrom,rs,gpos,ppos]
                rslist.append((chrom,int(ppos),rs))
        else:
             notfound += 1
             #print '###cannot find %s in ators' % aid
    print '# found',len(affyToMap),'affy ids and',notfound,'not found'
    rslist.sort()
    maplist = [(x[0],x[2],'0','%d' % x[1]) for x in rslist] # rslist was chrom, gpos,rs
    tabmaplist = ['\t'.join(x) for x in maplist]
    for race in ['white','black']:
        mapfname = os.path.join(gdir, '%s_%s.map' % (outfpref,race))
        mapf = file(mapfname,'w')
        mapf.write('\n'.join(tabmaplist))
        mapf.write('\n')
        print 'wrote %s' % mapfname
        mapf.close()
    return rslist,maplist


def getMarkers(annoroot='./', gdir='./', outfpref=None):
    """use affy's annotation files

[root@beast dbgap]# head -50 GenomeWideSNP_6.na23.annot.csv 
##For information about the Annotation file content, please see the bundled README file.
#%create_date=Wed Jul 18 10:52:51 PDT 2007
#%chip_type=GenomeWideSNP_6
#%chip_type=GenomeWideSNP_6.Full
#%genome-species=Homo_sapiens
#%genome-version=hg18
#%genome-version-ucsc=hg18
#%genome-version-ncbi=36.1
#%genome-version-create_date=2006-03
#%dbsnp-version=126
#%dbsnp-date=2006-05
#%netaffx-annotation-date=2007-07
#%netaffx-annotation-netaffx-build=23
#%netaffx-annotation-tabular-data-version=1.0
#%netaffx-annotation-tabular-format-version=1.0
"Probe Set ID","Affy SNP ID","dbSNP RS ID","Chromosome","Physical Position","Strand","ChrX pseudo-autosomal region 1","Cytoband","Flank","Allele A","Allele 
B","Associated Gene","Genetic 
Map","Microsatellite","Fragment Length Start Stop","Allele Frequencies","Heterozygous Allele Frequencies","Number of individuals/Number of chromosomes","In 
Hapmap","Strand Versus 
dbSNP","Copy Number Variation","Probe Count","ChrX pseudo-autosomal region 2","In Final List","Minor Allele","Minor Allele Frequency"
"SNP_A-1780419","---","rs6576700","1","84647761","-","0","p31.1","GGATACATTTTATTGC[A/G]CTTGCAGAGTATTTTT","A","G","ENST00000370662 // intron // 0 // --- // --- // --- 
// 
Deoxyribonuclease-2-beta precursor (EC 3.1.22.1) (Deoxyribonuclease II beta) (DNase II beta) (DNase2-like acid DNase) (DNase II- /// ENST00000370665 // intron // 0 
// --- // --- // --- // 
Deoxyribonuclease-2-beta precursor (EC 3.1.22.1) (Deoxyribonuclease II beta) (DNase II beta) (DNase2-like acid DNase) (DNase II- /// ENST00000260542 // intron // 0 
// --- // --- // --- // 
Deoxyribonuclease-2-beta precursor (EC 3.1.22.1) (Deoxyribonuclease II beta) (DNase II beta) (DNase2-like acid DNase) (DNase II- /// ENST00000361540 // intron // 0 
// --- // --- // --- // 
Deoxyribonuclease-2-beta precursor (EC 3.1.22.1) (Deoxyribonuclease II beta) (DNase II beta) (DNase2-like acid DNase) (DNase II- /// NM_058248 // intron // 0 // 
Hs.129142 // DNASE2B // 
58511 // Homo sapiens deoxyribonuclease II beta (DNASE2B), transcript variant 2, mRNA. /// NM_021233 // intron // 0 // Hs.129142 // DNASE2B // 58511 // Homo sapiens 
deoxyribonuclease II 
beta (DNASE2B), transcript variant 1, mRNA.","108.402547694204 // D1S2889 // D1S2766 // --- // --- /// 116.765222864103 // AFMA081XE9 // AFMB320YF1 // D1S2889 // 
D1S2766 /// 
102.95735950781 // --- // --- // 586552 // 615456","D1S1746E // upstream // 140114 /// D1S2889 // downstream // 39379","798 // 84647137 // 84647934 /// 1336 // 
84647550 // 84648885","0.6 
// 0.4 // Han Chinese /// 0.55 // 0.45 // Caucasian /// 0.6222 // 0.3778 // Japanese /// 0.5167 // 0.4833 // Yoruban","0.5778 // Han Chinese /// 0.5333 // Caucasian 
/// 0.5333 // Japanese 
/// 0.5 // Yoruban","45.0 // Han Chinese /// 60.0 // Caucasian /// 45.0 // Japanese /// 60.0 // Yoruban","YES","reverse","---","6","0","YES","G // Han Chinese /// G 
// Caucasian /// G // 
Japanese /// G // Yoruban","0.4 // Han Chinese /// 0.45 // Caucasian /// 0.3778 // Japanese /// 0.4833 // Yoruban"
"SNP_A-1780418","---","rs17054099","5","156323558","-","0","q33.3","GGATACATTACCCAAA[C/T]GGTCACAGGTCAAAGG","C","T","ENST00000274532 // upstream // 714 // --- // --- 
// --- // T-cell 
immunoglobulin and mucin domain-containing protein 4 precursor (TIMD-4) (T-cell membrane protein 4) (TIM-4). [Source:Unip /// ENST00000339252 // downstream // 65551 
// --- // --- // --- 
// Hepatitis A virus cellular receptor 1 precursor (HAVcr-1) (T cell immunoglobulin and mucin domain-containing protein 1) (TIMD-1) /// NM_138379 // upstream // 714 
// Hs.334907 // TIMD4 
// 91937 // Homo sapiens T-cell immunoglobulin and mucin domain containing 4 (TIMD4), mRNA. /// NM_012206 // downstream // 65551 // Hs.129711 // HAVCR1 // 26762 // 
Homo sapiens hepatitis 
A virus cellular receptor 1 (HAVCR1), mRNA.","161.231323294181 // D5S820 // D5S2049 // --- // --- /// 159.37910567766 // GATA51G10 // UT8057 // D5S1499 // D5S1403 
/// 154.762533280549 // 
--- // --- // 157925 // 272584","D5S2112 // upstream // 291875 /// D5S2417 // downstream // 196993","822 // 156323187 // 156324008","0.9667 // 0.0333 // Han Chinese 
/// 1.0 // 0.0 // 
Caucasian /// 1.0 // 0.0 // Japanese /// 0.825 // 0.175 // Yoruban","0.0667 // Han Chinese /// 0.0 // Caucasian /// 0.0 // Japanese /// 0.3167 // Yoruban","45.0 // 
Han Chinese /// 60.0 // 
Caucasian /// 45.0 // Japanese /// 60.0 // Yoruban","YES","reverse","---","6","0","YES","T // Han Chinese /// T // Caucasian /// T // Japanese /// T // 
Yoruban","0.0333 // Han Chinese /// 
0.0 // Caucasian /// 0.0 // Japanese /// 0.175 // Yoruban"
"
    """
    mapfname = os.path.join(gdir, '%s.map' % outfpref)
    try:
        mapf = file(mapfname,'r') # use existing map file - takes too long for testing otherwise
        ismap = 1
    except:
        ismap = 0
    if ismap:
        maplist = mapf.readlines()
        maplist = [x.strip().split() for x in maplist]
        rslist = [[x[0],x[3],x[1]] for x in maplist if len(x) > 1]
        print 'Read map',mapfname
        mapf.close()
    else:
       rslist = []
       rsSeen = {}
       fname = os.path.join(annoroot,'GenomeWideSNP_6.na23.annot.csv')
       crs= 'dbSNP RS ID' # csv reader column for rs number in affy annotation file
       cstart= 'Physical Position' # pompous pratt whoever named these columns
       cchrom= 'Chromosome'
       affyf = file(fname,'rb')
       for i in range(15): # yes, there are 15 useless lines
           crap = affyf.readline() # why oh why?
       f = csv.DictReader(affyf) # whew
       for n,row in enumerate(f):
            if n % 100000 == 0:
                print 'at row %d in %s' % (n,fname)
            rs = row[crs]
            if rs[0] <> '-':
                isold = rsSeen.get(rs,None)
                if isold:
                    print '### rs %s seen again - ignored' % rs
                    continue
                rsSeen.setdefault(rs,1)
                chrom=row[cchrom]
                start=row[cstart]
                try:
                    istart = int(start)
                    rslist.append((chrom,start,rs))
                except:
                    print 'rs %s start = %s' % (start,rs)
       rslist.sort()
       mapf = file(mapfname,'w')
       maplist = [(x[0],x[2],'0',x[1]) for x in rslist]
       tabmaplist = ['\t'.join(x) for x in maplist]
       mapf.write('\n'.join(tabmaplist))
       mapf.write('\n')
       print 'wrote %d rows to %s' % (len(tabmaplist),mapfname)
       mapf.close()
    return rslist,maplist


def getGenos(gf=None,rsdict=None,rsmissing={},nalleles=None):
    """read a subjects genos into the right order
       This version for dbGaP GAIN individual format
       Ross Lazarus Feb 2009
       about 25 secs on beast to ungzip and process one individual format file to a ped row

#Sample_ID: 11365
#Consent Group: General Research Use
ASSAY_SNP_ID ss_ID rs_ID strand orig_assay_orient_genotype genomic_fwd_orient_genotype confidence allele1_intensity allele2_intensity
SNP_A-2131660 ss76179473 rs2887286 + CC CC 5230 2913 573
SNP_A-1967418 ss76060847 rs1496555 + GG GG 4623 267 483
SNP_A-1969580 ss76062510 rs_undefined + GG GG 4325 725 3527
SNP_A-4263484 ss76379966 rs3890745 + TT TT 3800 313 1449
SNP_A-1978185 ss76068812 rs10492936 - CC GG 3713 1290 350
SNP_A-4264431 ss76380899 rs10489588 - CC GG 4161 1175 303
SNP_A-1980898 ss76070878 rs2376495 - GG CC 3515 365 837
    """
    adict = {'A':'1','C':'2','G':'3','T':'4','0':'0','N':'0','1':'1','2':'2','3':'3','4':'4'}
    # note Y appears in pseudo-autosomal X regions I think
    # not sure how to cope with - set to 0 for now - plink gets upset with that...

    grs = array.array('c','0'*nalleles) # all missing!
    for nrow,row in enumerate(gf):
        if row[0] <> '#':
            if len(row.strip().split()) > 5:
                lrow = row.strip().split()
                rs = lrow[2]
                rspos = rsdict.get(rs,None)
                if rspos and rspos < nalleles: # is a real rs
                    gt = lrow[5]
                    g = map(None,gt)
                    qc = lrow[6]
                    if (g <> ['0','0']) and lrow[-1] <> 'ND': # from share data - no idea if GAIN has this - leave in place
                        if g[1] == 'Y': # dbgap idea of a pseudo autosomal?
                           g[1] = g[0]  # hack to what we normally do for X         
                        # g = [adict.get(x,'0') for x in g] # translate into numbers ..
                        if g[0] == '0' or g[1] == '0': # otherwise plink barfs
                            g = ['0','0'] # wierd x psuedo-autosomal het crap?
                    grs[rspos] = g[0]
                    grs[rspos+1] = g[1]
                else:
                    if not rsmissing.get(rs,None):
                        rsmissing[rs] = rs
                        nm = len(rsmissing)
                        if debug:
                            if not rspos:
                                print '### rs for row=%s awol nm=%d' % (row,nm)
                            else:
                                print '### dud rspos %d > len(rsdict)*2 = %d (nalleles=%d) for row=%s' % (rspos,len(rsdict)*2,nalleles,row)
            else:
                print '## dud row = %s at %d' % (row,nrow)
    return grs,rsmissing
            
def readGenos(gdir = '/usr/dbgap/GAIN/allindfmt',outfpref='genos',genoprefix='*.geno.*'):
    """
    """
    controls,ped = readControls()
    outroot = os.path.join(gdir,outfpref)
    rslist,maplist = makeAffyMap(outfpref=outfpref,gdir=gdir)
    rss = [x[2] for x in rslist] # list of rs numbers
    nalleles = 2*len(rss)
    print 'nalleles = %d, len(rslist) = %d' % (nalleles,len(rslist))
    rssoffsets = range(0,nalleles,2)# index into array of alleles
    rsdict = dict(zip(rss,rssoffsets)) 
    sids = {}
    lrs = [x[2] for x in rslist] # just the rs numbers
    whiteoutf = file('%s_white.ped' % outroot,'w')
    blackoutf = file('%s_black.ped' % outroot,'w')
    bout = 0
    wout = 0
    rsmissing = {} # so only reported once
    started = time.time()
    for cnum,control in enumerate(controls.keys()):
        flist = glob.glob(os.path.join(gdir,'*%s.indfmt.*' % control))
        # note this avoids eg 61357.dup.indfmt.gz files
        if len(flist) == 0:
            print '### cannot find a file for control id %s' % control
            continue
        fname = flist[0]
        if os.path.splitext(fname)[-1] <> '.gz':
            print '### fname %s not a gzip' % fname
            continue
        p = ped.get(control,None)
        if p:
            c = controls.get(control,None)
            white = c['white'] == '1'
            gf = gzip.GzipFile(fname, 'rb');
            if ((cnum+1) % 10) == 0:
                dur = time.time() - started
                rate = dur/(cnum+1)
                print '%d secs: reading #%d = %s, %g secs each, %d white, %d black' % (dur,cnum+1,control,rate,wout,bout)
                blackoutf.flush()
                whiteoutf.flush()
            genos,rsmissing = getGenos(gf=gf,rsdict=rsdict,rsmissing=rsmissing,
            nalleles=nalleles)
            f = blackoutf
            if white:
                f = whiteoutf
                wout += 1
            else:
                bout += 1
            f.write(' '.join(p))
            f.write(' ')
            f.write(' '.join(genos))
            f.write('\n')
            gf.close()
        else:
            print '### control %s not found in pedigree data' % control
    blackoutf.close()
    whiteoutf.close()
    missout = 'rs_missing_from_JohnAnno.xls'
    print "# now writing all those mysterious rs numbers missing from John's Affy mapping to %s" % missout
    f = file(missout,'w')
    f.write('\n'.join(rsmissing.keys()))
    f.close()

if __name__ == "__main__":
    testMap = False
    if len(sys.argv) > 1:
        gdir = sys.argv[1].strip()
    else:
        gdir = '/usr/dbgap/GAIN/allindfmt'
    if len(sys.argv) > 2:
        outfpref = sys.argv[2]
    else:
        outfpref = 'GAINcontrolsx'
    if testMap:
        rslist,maplist = makeAffyMap(outfpref=outfpref,gdir=gdir) # jz fixed; ross refactored
        print 'rslist[100]=',rslist[:100]
        print 'maplist[100] = ',maplist[:100]
        raise 'OK'
    readGenos(gdir = '/usr/dbgap/GAIN/allindfmt',outfpref=outfpref,genoprefix='*.geno.*')


