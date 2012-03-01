"""
Affy intensities from framingham share individual genotype files
added January 2008 ross lazarus

Grrr - all records in individual format with "ND" for the various numeric values has
the first 4 fields gratuitously repeated so the parsing breaks unless the error is
patched.

note that you will need to edit /etc/security/limits.conf to raise the number of
allowed file handles (nofiles) to the number of subjects you want to process for the user who
will run this script - add lines for both hard and soft or it won't work. Limit is
changed at next login

getting the sort right is a major pain
this works for the original file
sort  -k2,2 -S8G -o Weiss_Harvard-Medical_FCR-DNA_Genotypes.rs Weiss_Harvard-Medical_FCR-DNA_Genotypes.sorted

sigh. 27k/sec is too slow. take 3 weeks for 600M rows.
grep for an individual rs takes too long - about 13 minutes

simply reading the file to find the "tell" values for each new rs
will take about 0.5 hours..

maybe just write files.
updated november 21 - will take weeks using current sort order
resorted into rs order so we can store one record per
marker with all the x,y and z values to plot

Weiss_Harvard-Medical_FCR-DNA_Genotypes.sorted header is

Sample ID       SNP Name        Allele1 - Forward       Allele2 - Forward       GC Score
Theta   R       X       Y       X Raw   Y Raw   B Allele FreLog R Ratio

Need to plot all data for each snp so must store snp wise - sqlite to the rescue?
will be a 60GB file - no worse I guess
will take a long time to make
"""
import time, sys, os, copy, subprocess, csv, glob

from sqlite3 import dbapi2 as sqlite


pch1 = '1' # arbitrary R plot characters
pch2 = '2'
pch3 = '3'
pch4 = '4'
swapalleles = {pch1:pch3,pch2:pch2,pch3:pch1,pch4:pch4} # swap pch dict
colours = {pch1:"dark blue", pch2:"dark green", pch3:"dark red", pch4:"salmon"}
plotsym = {pch1:1, pch2:3, pch3:2, pch4:4}
#colours = {pch1:"dark blue", pch2:"dark green", pch3:"dark red", pch4:"gray"}
debug = 0
mogresize = 'x300' # controls jpg thumbnail size


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

def timenow():    
   """return current time as a string    """    
   return time.strftime('%d/%m/%Y %H:%M:%S', time.localtime(time.time()))


    
class affyFCR:
    """Grab what we need from an affy 500k framingham individual genotype data report
    and stow in sqlite table for later recall and plotting
    unfortunately, the raw data is by subject - fortunately they're all in the same order
    so we read one line from each file and munge it into a single sqlite record for that rs
    
    Unfortunately, rs41388745 has two affy ids
    [rerla@hg affy500raw]$ grep rs41388745 *.13087
    ss66457619,rs41388745,+,na,AA,0,0,0.138336,1885.518066,651.891663
    ss66457618,rs41388745,+,na,AA,0,0,0.066873,2259.364014,709.739136
    
    """

    def __init__(self,dbname='/nfs/framinghamShare_nov2007/geno/affy500raw/affy500Framingham.sqlite',
                 tablename='intensities', gdir = '/nfs/framinghamShare_nov2007/geno/affy500raw',
                 genoprefix='*.geno.*'):
        """
        """
        self.gdir = gdir
        self.genoprefix = genoprefix
        self.dbname = dbname
        self.tablename = tablename
        self.intensityfile = os.path.split(dbname)[1] # filename
        self.varnames = ['id','rs','chrom','offset','a1','a2', 'GC','i1', 'i2', 'pch'] # 
        fieldnames = self.varnames # freshly read
        nfields = len(fieldnames)
        fieldtypes = ['text']*(nfields-1) # leave room for pk
        fieldtypes.insert(0,'INTEGER PRIMARY KEY NOT NULL') # make an autoincrement if
        self.fieldnames = fieldnames
        self.fieldtypes = fieldtypes
        self.createvars = ','.join(['%s %s' % (self.varnames[x],self.fieldtypes[x]) for x in range(nfields)])
        varplaces = ','.join(['?']*nfields) # sqlite uses qmarks
        self.insertsql = 'INSERT INTO %s (%s) values (%s)' % (self.tablename,','.join(self.varnames),varplaces)        
        self.con = sqlite.connect(self.dbname)
        self.cur = self.con.cursor()
 
    def getMarkers(self):
        """use affy's annotation files
    [rerla@meme affy500raw]$ head Mapping250K_Sty.na24.annot.csv
    ##For information about the Annotation file content, please see the bundled README file.
    "Probe Set ID","Affy SNP ID","dbSNP RS ID","Chromosome","Physical Position","Strand","ChrX pseudo-autosomal region 1","Cytoband","Flank","Allele 
    A","Al"
    "SNP_A-1780358","10001174","rs17325399","5","54442707","+","0","q11.2","ggataatacatattca[A/G]accaataacatacagc","A","G","NM_152623 // downstream 
    // 1873"
    "SNP_A-1780551","10003881","rs12454921","18","52924726","-","0","q21.31","ggatatcaaagttcca[A/G]gatgtcttttgaggtt","A","G","NM_052834 // 
    downstream // 76"


        """
        rslist = []
        fskel = '/nfs/framinghamShare_nov2007/geno/affy500raw/Mapping250K_%s.na24.annot.csv'
        crs= 'dbSNP RS ID' # csv reader column for rs number in affy annotation file
        cstart= 'Physical Position' # pompous pratt whoever named these columns
        cchrom= 'Chromosome'
        for chip in ['Sty','Nsp']:
            fname = fskel % chip
            affyf = file(fname,'rb')
            crap = affyf.readline() # why oh why?
            f = csv.DictReader(affyf) # whew
            for n,row in enumerate(f):
                if n % 100000 == 0:
                    print 'at row %d in %s' % (n,fname)
                rs = row[crs]
                chrom=row[cchrom]
                start=row[cstart]
                try:
                    start = int(start)
                    rslist.append((chrom,start,rs))
                except:
                    pass #print 'rs %s start = %s' % (start,rs)
        rslist.sort()
        return rslist
            

    def readPed(self,pedname="/nfs/framinghamShare_nov2007/geno/affy500raw/phs000007.pht000183.v1.p1.c2.share_ped.NPU.txt",
                ped2name='/nfs/framinghamShare_nov2007/geno/affy500raw/phs000007.pht000182.v1.p1.c2.shareids.NPU.txt'):
        """
        make the start of a pedigree file
    # Study accession: phs000007.v1.p1
    # Table accession: pht000183.v1.p1.c2
    # Consent group: Non-Profit Use Only
    # Citation instructions: The study accession (phs000007.v1.p1) is used to cite the study and its data tables and documents. The data in this 
    file shoul 
    cited us$
    # To cite columns of data within this file, please use the variable (phv#) accessions below:
    #
    # 1) the table name and the variable (phv#) accessions below; or
    # 2) you may cite a variable as phv#.v1.p1.c2.

    ##phv00024067   phv00024068     phv00024069     phv00024070     phv00024071     phv00024072     phv00024073     phv00024074     phv00024075     
    phv0002 
    phv0$
    pedno   shareid fshare  mshare  SEX     itwin   idtype  pheno   geno    genrelease      nonprofit
    1       16                      2               0       1               0       1
    1       3226    9218    8228    2               3       1       1       1       0
    1       8228                    2               1       1       1       1       0
    1       9218    17008   16      1               1       1               1       0


        not families, but this file contains lots of subjects missing from the share_ped file!
        [rerla@meme geno]$ more phs000007.pht000182.v1.p1.c2.shareids.NPU.txt 
    # Study accession: phs000007.v1.p1
    # Table accession: pht000182.v1.p1.c2
    # Consent group: Non-Profit Use Only
    # Citation instructions: The study accession (phs000007.v1.p1) is used to cite the study and its data tables and documents. The data in this 
    file should
     be cited using the accession pht000182.v1.p1.c2.
    # To cite columns of data within this file, please use the variable (phv#) accessions below:
    #
    # 1) the table name and the variable (phv#) accessions below; or
    # 2) you may cite a variable as phv#.v1.p1.c2.

    ##phv00024058   phv00024059     phv00024060     phv00024061     phv00024062     phv00024063     phv00024064     phv00024065     phv00024066
    shareid idtype  SEX     pheno   geno    pedno   itwin   genrelease      nonprofit
    1       0       2       1                               0       1
    2       1       2       1               126             0       1
    3       0       2       1       1                       0       1

        """
        f = file(pedname,'r')
        header = None
        ped = {}
        for row,line in enumerate(f):
            line = line.strip()
            if len(line) > 0:
                if line[0] <> '#':
                    if header == None:
                        header = line.split('\t')
                        cfamid = header.index('pedno')
                        cshareid = header.index('shareid')
                        cfshareid = header.index('fshare')
                        cmshareid = header.index('mshare')
                        cgender = header.index('SEX')
                    else:
                        ll = line.split('\t')
                        pednum = ll[cfamid]
                        shareid = ll[cshareid]
                        fid = '0'
                        if len(ll[cfshareid]) > 0:
                            fid = ll[cfshareid]
                        mid = '0'
                        if len(ll[cmshareid]) > 0:
                            mid = ll[cmshareid]
                        gender = ll[cgender]
                        ped[shareid] = [pednum,shareid,fid,mid,gender,'1'] # dummy affection status
        print 'read %d rows from %s' % (row,pedname)
        f = file(ped2name,'r')
        header = None
        ped = {}
        for row,line in enumerate(f):
            line = line.strip()
            if len(line) > 0:
                if line[0] <> '#':
                    if header == None:
                        header = line.split('\t')
                        cfamid = header.index('pedno')
                        cshareid = header.index('shareid')
                        cgender = header.index('SEX')
                    else:
                        ll = line.split('\t')
                        shareid = ll[cshareid]
                        pednum = ll[cfamid]
                        if len(pednum) == 0:
                            pednum = '999%s' % shareid
                        fid = '0'
                        mid = '0'
                        gender = ll[cgender]
                        ped[shareid] = [pednum,shareid,fid,mid,gender,'1'] # dummy affection status
        print '##read %d rows from %s' % (row,pedname)
        print '##found %d pedigrees' % (len(ped.keys()))
        return ped

    def getMarkerdat(self,rsnum='rs1234'):
        """ return all data for this marker in a useful way
        """
        sql = 'select * from %s where rs = "%s"' % (self.tablename,rsnum)
        self.cur.execute(sql)
        res = self.cur.fetchall()
        if (len(res) > 0):
            return res[0]
        else:
            return None

    def getRegion(self,chrom='1',spos=100000,epos=1000000):
        """lookup region
        """
        sql = 'select rs from %s where chrom = "%s" and round(offset) > %d and round(offset) <= %d' % (self.tablename,
                                                                        chrom,spos,epos)
        self.cur.execute(sql)
        res = self.cur.fetchall()
        if len(res) > 0:
            res = [x[0] for x in res] # list of rs
        else:
            res = []
        return res
            
    def intensityIter(self):
        """Affy version for framingham individual file
        format - only one with intensity values.
        Take a line from each file and create the kind of record
        mysql and sqlite expect for executemany
        Inserts a Null value for primary key which sqlite will autoincrement
        grs.append(' '.join((rs,id,gt,qc,a1i,a2i))) # for quickly writing out
        
        rerla@meme affy500raw]$ head phg000006.ind.geno.NonProfit_Release.geno.10481
        #File_Name: phg000006.ind.geno.NonProfit_Release.geno.10481
        #Consent_group: Non-profit_use
        #Sample_ID1: 10481
        #ss#,rs#,ss2rsOrientation,rs2GenomeOrientation,genotype_genomic_orient,genotye_orig_code,
        qc_bit_flags,genotype_quality_score,allele1_intensity,allele2_y
        ss66079302,rs3094315,+,-,GA,1,0,0.002786,482.392090,445.126862
        ss66273559,rs4040617,+,+,AG,1,0,0.006784,768.500610,561.713074
        ss66317030,rs2980300,-,+,TC,1,0,0.001729,388.633301,387.334625
        ss66185183,rs2905036,-,-,TT,2,0,0.031937,196.904343,703.350342
        ss66174584,rs4245756,+,+,CC,0,0,0.018019,1021.935425,206.439590
        ss66145570,rs4075116,+,-,TC,1,0,0.002691,888.154358,673.292847
     
        """
        # chrom,start,rs
        listrs = [x[2] for x in self.rslist]
        self.rsdict = dict(zip(listrs,self.rslist)) # chrom, start, rs keyed by rs
        self.ped = self.readPed()
        gglob = os.path.join(self.gdir,self.genoprefix)
        gflist = glob.glob(gglob) # should be about 9200 files or so
        print '## %s -> len gflist=%d' % (gglob,len(gflist))
        outroot = 'framinghamaffy'
        listf = []
        ids = []
        for i,fname in enumerate(gflist): # check that we have a pedigree for each subject
            shareid = fname.split('.')[-1] # should be the very last bit
            p = self.ped.get(shareid,None)
            if p:
                try:
                    f = file(fname,'r')
                except:
                    print '### Mmmm. Problem opening file # %d called %s' % (i,fname)
                    print '### Do you need to increase the number of open files (nopenf)'
                    print '### in /etc/security/limits.conf?'
                    print '### or is this a file read permission error?'
                    sys.exit(1)
                listf.append(f) # make a list of file objects - one per subject
                ids.append(shareid)
            else:
                print '## ped has no shareid %d' % shareid
        print '###listf has %d files' % (len(listf))
        self.listf = listf
        batchsize = 1000
        wewant = ['rs','gt', 'gqc','i1','i2'] # this determines which fields we keep!
        wewantcols = [1,4,7,8,9]
        lnum = 0
        currentrs = None # sort order is rs so we can batch
        res = []
        done = 0
        writtenrs = {}
        nrecs = 0
        started = time.time()
        while not done:
            rsa1 = rsa2 = None # reinitialise all the accumulators for this new marker
            rsx = []
            rsy = []
            rsgc = []
            rspch = []
            alleles = {}
            firsta1 = None
            currentrs = None
            lnum += 1
            for fid,ifile in enumerate(self.listf): # each time we cycle through, we have all data for a marker
                try:
                    line = ifile.readline()
                    while line[0] == '#':
                        line = ifile.readline()
                except:
                    done = 1
                    print 'done'
                    if len(res) > 0:
                        yield res
                    raise StopIteration
                ll = line.strip().split(',')
                if len(ll) > 10:
                    ll = ll[4:] # ND records are b0rken in the current data = first 4 fields repeated
                rs,gt,qc,x,y = [ll[x] for x in wewantcols] # filter only the fields we want to save
                if (rs <> currentrs): # infiles MUST be sorted by rs for this to work
                    if currentrs == None: # first time
                        currentrs = rs
                        chrom,offset,rsa = self.rsdict.get(rs,('?',0,'rs0'))
                        if rsa <> rs:
                            print '### rsdict lookup of %s got %s' % (rs,rsa)
                        if writtenrs.get(currentrs,None):
                            # problem - we wrote this rs already
                            print '### lnum%d: rs=%s, currentrs=%s and already written...' % (lnum,rs,currentrs)
                        writtenrs[currentrs] = currentrs
                    else: # we have a problem!
                        print 'oh dear, line %d, file# %d has rs %s instead of rs %s' % (lnum,i,rs,currentrs)
                        sys.exit(1)
                a1 = gt[0]
                a2 = gt[1]
                if a1 == 'N' or a2 == 'N' or ll[-1] == 'ND': # missing = b0rked records!
                    plotchar = pch4 # plot a missing
                elif a1 == 's' or a2 == 's': # mystery
                    plotchar = pch4
                    print '##WTF: for id %s, we got ll %s' % (ids[fid],ll)
                else:
                    alleles[a1] = alleles.get(a1,0) + 1
                    alleles[a2] = alleles.get(a2,0) + 1
                    if firsta1 == None:
                        firsta1 = a1
                    if a1 == a2: # hom
                        if a1 == firsta1:
                            plotchar = pch1 # plot as common hom.
                        else:
                            plotchar = pch3 # plot as rare hom.
                    else:
                        plotchar = pch2
                rspch.append(plotchar)
                try:
                    test = float(x)
                    rsx.append(x) # make vectors in subject order - don't really care who
                except:
                    rsx.append('0.0')
                try:
                    test = float(y)
                    rsy.append(y)
                except:
                    rsy.append('0.0')
                try:
                    test = float(qc)
                    rsgc.append(qc)
                except:
                    rsgc.append('0.0')
            # at end of subject files list for one rs
            if len(rspch) > 0:
                if len(alleles.keys()) <> 2:
                    print '## ignoring %s because we see alleles = %s' % (currentrs,alleles)
                    # this happens eg for monomorphs rs1288356 that we're going to ignore...
                else:
                    
                    lallele = alleles.keys() # the alleles
                    rsa1,rsa2 = lallele[:2]
                    if alleles[rsa2] > alleles[rsa1]: # need to swap so rsa1 is common
                        rsa1,rsa2 = rsa2,rsa1
                    if rsa1 <> firsta1:
                        rspch = [swapalleles[x] for x in rspch] # consistent order for major allele
                    if (len(alleles.keys()) <> 2) or (rsa1=='s') or (rsa2=='s'):
                        print 'for rs %s we have alleles = %s, rsa1=%s,rsa2=%s' % (currentrs,alleles,rsa1,rsa2)
                    res.append([None,currentrs,chrom,offset,rsa1,rsa2,','.join(rsgc),','.join(rsx),','.join(rsy),','.join(rspch)])
                    if len(res) >= batchsize:
                        nrecs += len(res)
                        dur = time.time() - started
                        rate = nrecs/dur
                        print '###yielded %d records @ %f per second' % (nrecs,rate)
                        yield res
                        res = []
            else:
                print '### empty rspch'
                done = 1
        if len(res) > 0: # last few records if not already caught
            yield res

            
    def buildSQLite(self):
        """create and fill a sqlite db from the mapping annotation files downloaded from
        http://www.affymetrix.com/support/technical/byproduct.affx?product=500k
        """
        # Create a connection to the database file "mydb"
        self.rslist = self.getMarkers()
        try:
            sql = 'drop table %s' % self.tablename
            self.cur.execute(sql)
        except:
            if debug:
                print 'no %s to drop this time!' % self.tablename
            else:
                pass
        sql = 'create table %s (%s)' % (self.tablename,self.createvars)
        print 'trying %s' % sql
        self.cur.execute(sql)
        print 'executed %s' % sql
        sql = 'create index rsi on %s (rs)' % self.tablename
        self.cur.execute(sql)
        dat = self.intensityIter()
        for recs in dat: # generator is iterable - ah, python.
            if len(recs) == 0:
                print '## recs empty in buildSQLite'
                break
            self.cur.executemany(self.insertsql,recs) #
            self.con.commit()
        self.con.commit() # final cleanup
        
    def buildMySQL(self,host, userid, password):
        """create and fill a mysql db from the mapping annotation files downloaded from
        http://www.affymetrix.com/support/technical/byproduct.affx?product=500k
        """
        import MySQLdb
        ntoinsert = 500
        # Create a connection to the database file "mydb"
        genome = MySQLdb.Connect('godzilla', 'refseq', 'Genbank')
        curs = genome.cursor() # use default cursor
        curs.execute('use %s' % self.dbname)
        varss = '''id INT UNSIGNED NOT NULL AUTO_INCREMENT, affyid varchar(15),probesetid varchar(15),
        rs varchar(15), chrom char(2), offset int(12), strand char(1),     
        index rsindex (rs), index chromrs (chrom, rs), primary key (id)'''
        ivar = "id,affyid,probesetid,rs,chrom,offset,strand"
        qvall = ['%s']*7
        qvals = ','.join(qvall)
        insertsql = """insert into %s (%s) values (%s)""" % (self.tablename,ivar,qvals)
        sql = 'DROP TABLE IF EXISTS %s.%s' % (self.dbname,self.tablename) # start again each time
        curs.execute(sql)
        sql = 'CREATE TABLE %s.%s (%s)' % (self.dbname,self.tablename,varss)      
        curs.execute(sql)
        f = affyannoF(self.infname)
        done = 0
        insertme = []
        while not done:
                while (len(insertme) < ntoinsert) and not done:
                    try:
                        ll = f.next()
                        insertme.append(ll)
                    except StopIteration:
                        done = 1
                if len(insertme) > 0:
                    curs.executemany(insertsql,insertme) # just pass the list
                    insertme = []
                genome.commit() # only at end of each file to speed things up
        curs.close()
        genome.close()

def fail( message ):
    print >> sys.stderr, message
    return -1
   

def makeplots(rslist=['rs1234',],title='title',outfile='test',
              h=8,w=8,outdir='./', xyfname='xy',mogrify=0,region='',progname='foo',nup=4):
    """
    y is data for a plot of genotype call by x and y values
    from an affy full call report
    note can layout symmetrical figures as
    foo = layout(matrix(c(1:9),nrow=3,byrow=T)) but must specify divisible row or col
    so probably best to make separate outputs and put them in an html style galaxy output?
    """
    idb = affyFCR(dbname='/mnt/memefs/framinghamShare_nov2007/geno/affy500raw/affy500Framingham.sqlite',
                 tablename='intensities', gdir = '/mnt/memefs/framinghamShare_nov2007/geno/affy500raw',
                  genoprefix='*.geno.*')
    if region > '':
        try: # TODO make a regexp?
            c,rest = region.split(':')
            c = c.replace('chr','')
            try:
                c = int(c)
            except:
                pass
            rest = rest.replace(',','') # remove commas
            spos,epos = rest.split('-')
            spos = int(spos)
            epos = int(epos)
            print >> sys.stdout, '## %s parsing chrom %s from %d to %d' % (progname, c,spos,epos)
        except:
            print >> sys.stdout, '##! %s unable to parse region %s - MUST look like "chr8:10,000-100,000' % (progname, region)
            sys.exit(1)
        rslist = idb.getRegion(chrom=c,spos=spos,epos=epos)
        if rslist == None or len(rslist) == 0:
            print >> sys.stdout, '##! %s no snps found in region %s' % (progname, region)
            sys.exit(1)
        else:
            rslist = [x.encode('ascii') for x in rslist]
            print '### got rslist=',rslist
    else:
        rsdict = dict(zip(rslist,rslist)) # for quick lookups
        print >> sys.stdout, '##%s looking for %d rs (%s)' % (progname,len(rslist),rslist[:5])
    intensityfile=idb.intensityfile
    flist = [] # will become urls for pdfs in html file output
    explanations = [] # link titles for pdfs
    rslist.sort(lambda x,y: int(x[2:])-int(y[2:])) # numeric rs sort 
    import rpy
    for marker in rslist:
        minfo = idb.getMarkerdat(marker)
        if minfo <> None: # exists
            #res.append([None,currentrs,rsa1,rsa2,','.join(rsgc),','.join(rsx),','.join(rsy),','.join(rspch)])
            id,rs,chrom,offset,a1,a2,gcvec,xvec,yvec,pchvec = minfo # unpack
            chrom = chrom.encode('ascii')
            offset = offset.encode('ascii')
            l1 = '%s/%s' % (a1,a1)
            l2 = '%s/%s' % (a1,a2)
            l3 = '%s/%s' % (a2,a2)
            l4 = 'Missing'
            llist = [l1,l2,l3,l4]
            llist = [x.encode('ascii') for x in llist]
            plist = [pch1,pch2,pch3,pch4]
            clist = [colours[x] for x in plist]
            plist = [plotsym[x] for x in plist]
            xvec = [float(i) for i in xvec.split(',') if i <> 'ND']
            yvec = [float(i) for i in yvec.split(',') if i <> 'ND']
            zvec = [int(i) for i in pchvec.split(',')]
            colvec = [colours['%d' % i] for i in zvec]
            pchvec = [plotsym['%d' % i] for i in zvec]
            nrows = len(xvec)
            subt = 'Maj/Min=%s/%s (N=%d)' % (a1,a2,nrows)
            main = '%s. %s chr%s:%s' % (title,marker,chrom,offset)
            ylab = 'Allele 2 intensity'  
            xlab = 'Allele 1 intensity'	
            fname = '%s_%s.pdf' % (intensityfile,marker)
            flist.append(fname)
            expl = '%s Intensity Plot' % marker
            explanations.append(expl)
            rpy.r.pdf( os.path.join(outdir,fname), h , w  )
            rpy.r.par(mai=(1,1,1,0.5))
            rpy.r.plot(xvec,yvec,xlab=xlab,ylab=ylab,main=main, pch=pchvec, col=colvec)
            rpy.r.grid(nx = None, ny = None, col = "lightgray", lty = "dotted")
            rpy.r.legend("topright",legend=llist,pch=plist,col=clist,title="Call")
            rpy.r.dev_off()
        else:
            flist.append('')
            explanations.append('##ERROR: %s NOT found in %s database##' % (marker,intensityfile))
    if nup:
        # pdfjoin --outfile chr1test.pdf `ls database/files/dataset_396_files/*.pdf`
        # pdfnup chr1test.pdf --nup 3x3 --frame true --outfile chr1test3.pdf 
        filestojoin = ' '.join(flist) # all the file names so far
        afname = '%s_All_Paged.pdf' % (intensityfile)
        flist.insert(0,afname) # bump
        expl = 'All %s Intensity Plots joined into a single pdf' % intensityfile
        explanations.insert(0,expl)
        vcl = 'pdfjoin %s --outfile %s ' % (filestojoin, afname)
        # make single page pdf  
        x=subprocess.Popen(vcl,shell=True,cwd=outdir)
        retval = x.wait()
        nfname = '%s_All_%dx%d.pdf' % (intensityfile,nup,nup)
        flist.insert(0,nfname)
        expl = 'All %s Intensity Plots %d by %d to a page' % (intensityfile,nup,nup)
        explanations.insert(0,expl)
        vcl = 'pdfnup %s --nup %dx%d --frame true --outfile %s' % (afname,nup,nup,nfname)
        # make thumbnail images    
        x=subprocess.Popen(vcl,shell=True,cwd=outdir)
        retval = x.wait()
    if mogrify:
        vcl = 'mogrify -format jpg -resize %s %s' % ( mogresize, '*.pdf')
        # make thumbnail images    
        x=subprocess.Popen(vcl,shell=True,cwd=outdir)
        retval = x.wait()
        

    return flist,explanations
    
def plotMaker():
    """
    <command interpreter="python2.4">
        affyPlotmaker.py "$rslist" $outfile1 $outfile1.files_path "$title" $height $width $mogrify "$region" "$nup"
    </command>

        """
    progname = os.path.split(sys.argv[0])[-1] # get script name
    print >> sys.stdout,'## %s at %s cl=%s' % (progname, timenow(), ' '.join(sys.argv))
    orslist = sys.argv[1].lower().replace('__cr__',' ')
    orslist = orslist.replace('__cn__',' ').split() # galaxy replaces newlines - go figure
    badmarkers = []
    region = ''
    if len(sys.argv) > 8:
        region = sys.argv[8]
    if len(sys.argv[1]) > '' and region == '':
        rslist = []
        for marker in orslist: # check that we have legal rsN..N markers so sort will work
            bad = 0
            try:
                x = int(marker[2:])
            except:
                bad = 1
            if not bad:
                rslist.append(marker)
            else:
                badmarkers.append(marker)
    else:
        rslist = ''
    outfile = sys.argv[2]
    newfilepath = sys.argv[3]
    title = sys.argv[4] # unchanged
    height = int(sys.argv[5])
    width = int(sys.argv[6])
    thumbnails = (sys.argv[7].lower() == "true")
    nup = None
    if len(sys.argv) > 9:
        try:
            nup = int(sys.argv[9])
        except:
            pass
    xyfname = 'X_Y'
    try:
        os.makedirs(newfilepath)
    except:
        pass    
    newfiles,explanations = makeplots(rslist=rslist,title=title,outfile=outfile,outdir=newfilepath,
              xyfname=xyfname,h=height,w=width,mogrify=thumbnails,region=region,progname=progname,nup=nup)
    for i in range(len(badmarkers)):
        newfiles.insert(0,'')
        explanations.insert(0,'##ERROR: Illegal marker %s! "rs" followed by digits ONLY##' % badmarkers[i])
    lf = file(outfile,'w')
    lf.write(galhtmlprefix % progname)
    s = '<div>Output from %s run at %s<br>\n' % (progname,timenow())       
    lf.write('<h4>%s</h4>\n' % s)
    if thumbnails: # tabular output
        lf.write('<div><h4>Click the images below to view PDF Intensity plots</h4>\n')        
        lf.write('<table cellpadding="5" border="0">')
    else:
        lf.write('<div><h4>Click the links below to view Intensity plots for each marker</h4><br><ol>\n')
    for i in range(len(newfiles)):
        title = explanations[i]
        url = newfiles[i]
        if url > '':
            if thumbnails:
                imghref = os.path.splitext(url)[0] # removes .pdf
                lf.write('<tr><td><a href="%s"><img src="%s.jpg" alt="%s" hspace="10" align="middle"></a>' \
                    % (url,imghref,title))
                lf.write('</td><td><a href="%s">%s</a></td></tr>\n' % (url,title))
            else:    
                lf.write('<li><a href="%s">%s</a></li>\n' % (url,title))
        else:
            if thumbnails:
                lf.write('<tr><td>&nbsp;</td><td>%s</td></tr>' % (title))
            else:
                lf.write('<li>%s</li>\n' % (title))
    if thumbnails:
        lf.write('</table><hr>\n')
    else:
        lf.write('</ol></div>\n')
    lf.write('</div></body></html>\n')
    lf.close()
    print >> sys.stdout, '%s closed' % outfile



if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1].lower() == 'buildaffy':
        foo = affyFCR(dbname='/mnt/memefs/framinghamShare_nov2007/geno/affy500raw/affy500Framingham.sqlite',
                 tablename='intensities',gdir = '/mnt/memefs/framinghamShare_nov2007/geno/affy500raw',
                    genoprefix='*.geno.*')
 
        foo.buildSQLite() # be warned, this takes hours!
    else:
        print >> sys.stdout, "cl=%s" % (' '.join(sys.argv))
        plotMaker()

