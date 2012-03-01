"""
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

feb 21, ross: ah. gc is a property of each assay, not each subject. removed redundant array.

"""
import time, sys, os, copy, Numeric, subprocess, csv, glob

from pysqlite2 import dbapi2 as sqlite

fileheader = ['Sample ID','SNP Name','Allele1 - Forward','Allele2 - Forward','GC Score',
'Theta','R','X','Y','X Raw','Y Raw','B Allele Freq','Log R Ratio']
pch1 = '1' # arbitrary R plot characters
pch2 = '2'
pch3 = '3'
pch4 = '4'
swapalleles = {pch1:pch3,pch2:pch2,pch3:pch1,pch4:pch4} # swap pch dict
colours = {pch1:"dark blue", pch2:"dark green", pch3:"dark red", pch4:"black"}
plotsym = {pch1:1, pch2:3, pch3:2, pch4:4}
#colours = {pch1:"dark blue", pch2:"dark green", pch3:"dark red", pch4:"gray"}
debug = 1
mogresize = 'x200' # controls jpg thumbnail size


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

class illFCR:
    """Grab what we need from an illumina infinium full call report
    and stow in sqlite table for later recall and plotting
    """

    def __init__(self,fpath='/usr/local/galaxy/data/camp2007',dbname='illcamp.sqlite',tablename='intensities',
     infile='Weiss_Harvard-Medical_FCR-DNA_Genotypes.rs'):
        """
        """
        self.dbname = os.path.join(fpath,dbname)
        self.tablename = tablename
        self.infile = os.path.join(fpath,infile)
        self.varnames = ['id','rs','chrom','offset','a1','a2', 'GC','Theta','R','X','Y',
                         'rawX','rawY','BAF','LogRR','pch'] # 
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
        if debug:
            print '++++ill initialised+++++'
            print 'got fieldnames = %s' % fieldnames
            print 'createvars = %s' % self.createvars
            print 'insertsql = %s' % self.insertsql
            print '++++ill initialised+++++'

    def getMap(self,mapf='/usr/local/galaxy/data/rg/library/lped/camp2007.map',
               annof='/home/rerla/camp2007/HumanHap550v3_A.csv'):
        """first priority for truth. Use original illumina anno as secondary source

        [rerla@hg ~]$ head camp2007/HumanHap550v3_A.csv 
        "Illumina, Inc. 2006"
        [Heading]
        Descriptor File Name,HumanHap550v3_A.csv
        Assay Format,Infinium II
        Date Manufactured,12/21/2006
        Loci Count,561466
        [Assay]
        IlmnID,Name,IlmnStrand,SNP,AddressA_ID,AlleleA_ProbeSeq,AddressB_ID,AlleleB_ProbeSeq,GenomeBuild,Chr,MapInfo,Ploidy,Species,Source,SourceVersion,Sourcd
        MitoA10045G-13273284_B_R_IFB1141652022:0,MitoA10045G,Bot,[T/C],0904860368,GAGGGTGTTGATTATTAAAATTAAGGCGAAGTTTATTACTCTTTTTTGAA,,,AF347015.1,chrM,10045,11
        MitoA10551G-13273286_T_F_IFB1141639111:        
        """
        try:
            f = open(mapf,'r')
        except:
            raise '%s not found or not openable - permission or existential problem'
        map = {}
        for n,row in enumerate(f):
            row = row.strip().split()
            if len(row) == 4 and n > 0:
                chrom,rs,gmap,offset = row
                map[rs] = row
        print 'read %d snps from %s' % (n,mapf)

        wewant = ['Name','IlmnID','SNP','Chr','MapInfo'] # this determines which fields we keep!
        f = file(annof,'r')
        header = ''
        while header.strip() <> '[Assay]':
            header = f.next()
        header = f.next()
        fieldnames = header.strip().split(',') # is csv
        fieldnames = [x.replace(' ','_') for x in fieldnames]
        fieldnames = [x.replace('-','') for x in fieldnames]
        fieldnames = [x.replace('"','') for x in fieldnames] # get rid of quotes
        header = fieldnames
        wewantin = [fieldnames.index(x) for x in wewant] # keep these fields
        if debug:
            print 'header=%s' % fieldnames
            print 'wewant=%s' % wewant
        minrow = max(wewantin)
        for row in f:
            row = row.strip().split()
            if len(row) >= minrow:
                fields = [row[x] for x in wewantin] # mmmm
                rs,illid,alleles,chrom,offset = fields
                if not map.get(rs,None): # use illumina universe
                    map[rs] = [chrom,rs,0,offset]
        print 'getmap returning map with %d rs' % (len(map.keys()))                    
        return map
                
    
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
            
    def illiter(self):
        """ heavily modified to deal with rs sorted file
        and to make a single db entry for each rs
        iterable for tab file lines as field lists
        mysql and sqlite expect for executemany
        Inserts a Null value for primary key which sqlite will autoincrement
        
        """
        map = self.getMap()
        batchsize = 1000

        header = ['Sample ID','SNP Name','Allele1 - Forward','Allele2 - Forward','GC Score','Theta','R','X','Y',
                'X Raw','Y Raw','B Allele Freq','Log R Ratio']
        wewant = ['SNP Name','Allele1 - Forward','Allele2 - Forward',
                     'GC Score','Theta','R','X','Y','X Raw','Y Raw','B Allele Freq','Log R Ratio']
        # this determines which fields we keep!
        f = file(self.infile,'r')
        print 'header=',header
        print 'wewant=',wewant
        wewantin = [header.index(x) for x in wewant] # keep these fields
        print 'wewantin=%s' % wewantin
        minlen = max(wewantin)
        lnum = 1
        started = time.time()
        currentrs = None # sort order is rs so we can batch
        res = []
        done = 0
        writtenrs = {}
        for lnum, line in enumerate(f):
            if line[:6] == 'Sample':
                print '### eating more header...'
            else:
                ll = line.strip().split('\t')
                if lnum % 1000000 == 0:
                    dur = time.time() - started
                    self.con.commit()
                    print 'illiter reading line %d of %s = %5.2f recs/sec' % (lnum,self.infile,lnum/dur)
                if len(ll) >= minlen:
                    ll = [ll[x] for x in wewantin] # filter only the fields we want to save
                    rs,a1,a2,gc,theta,r,x,y,rx,ry,baf,lrr = ll
                    if (rs <> currentrs): # infile MUST be sorted by rs for this to work
                        if currentrs <> None: # not first time - append data to res
                            if writtenrs.get(currentrs,None):
                                # problem - we wrote this rs already
                                print '### lnum%d: rs=%s, currentrs=%s and already written...' % (lnum,rs,currentrs)
                            writtenrs[currentrs] = currentrs
                            acounts = alleles.values()
                            allele = alleles.keys() # the alleles
                            if len(allele) < 2: # ignore lots!
                                acounts += [0,0]
                                allele += ['-','-']
                            rsa1 = allele[0]
                            rsa2 = allele[1]
                            if acounts[1] > acounts[0]: # need to swap
                                rsa1 = allele[1]
                                rsa2 = allele[0]
                            if rsa1 <> firsta1:
                                rspch = [swapalleles[x] for x in rspch] # consistent order for major allele
                            m = map.get(currentrs,None)
                            if m:
                                chrom,mrs,goffset,offset = m
                                res.append([None,currentrs,chrom,offset,rsa1,rsa2,rsgc,','.join(rstheta),
                                       ','.join(rsr),','.join(rsx),','.join(rsy),','.join(rsrx),','.join(rsry),
                                            ','.join(rsbaf),','.join(rslogrr),','.join(rspch)])
                            if len(res) >= batchsize:
                                yield res
                                res = []
                        rsa1 = rsa2 = None # reinitialise all the accumulators for this new marker
                        rstheta = []
                        rsx = []
                        rsy = []
                        rsr = []
                        rsrx = []
                        rsry = []
                        rsbaf = []
                        rslogrr = []
                        rspch = []
                        alleles = {}
                        firsta1 = None
                        rsgc = gc
                        currentrs = rs
                    if a1 == '-' or a2 == '-': # missing
                        plotchar = pch4 # plot a missing
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
                    rsx.append(x) # make vectors in subject order - don't really care who
                    rsy.append(y)
                    rsr.append(r)
                    rsrx.append(rx)
                    rsry.append(ry)
                    rsbaf.append(baf)
                    rslogrr.append(lrr)
                    rstheta.append(theta)
                else:
                    print 'line %d too short - %s' % (lnum,ll)
        if len(res) > 0:
            yield res
        raise StopIteration                

            
    def buildSQLite(self):
        """create and fill a sqlite db from the mapping annotation files downloaded from
        http://www.affymetrix.com/support/technical/byproduct.affx?product=500k
        """
        # Create a connection to the database file "mydb"
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
        dat = self.illiter()
        for recs in dat: # generator is iterable - ah, python.
            if len(recs) == 0:
                print '## recs empty in buildSQLite'
                break
            self.con.executemany(self.insertsql,recs) # 
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
   

def makeplots(mlist=['rs1234',],title='title',outfile='test',
              h=6,w=8,outdir='./', xyfname='xy',rthetafname='rtheta',mogrify=0):
    """
    y is data for a plot of genotype call by x and y values
    from an illumina full call report
    note can layout symmetrical figures as
    foo = layout(matrix(c(1:9),nrow=3,byrow=T)) but must specify divisible row or col
    so probably best to make separate outputs and put them in an html style galaxy output?
    """
    import rpy
    cex = 0.3
    idb = illFCR()
    #idb = illFile()
    flist = [] # will become urls for pdfs in html file output
    explanations = [] # link titles for pdfs
    mlist.sort(lambda x,y: int(x[2:])-int(y[2:])) # numeric rs sort 
    for marker in mlist:
        minfo = idb.getMarkerdat(marker)
        if minfo <> None: # exists
            id,rs,chrom,offset,a1,a2,gc,thetavec,rvec,xvec,yvec,rxvec,ryvec,rbaf,rlrr,pchvec = minfo # unpack
            l1 = '%s/%s' % (a1,a1)
            l2 = '%s/%s' % (a1,a2)
            l3 = '%s/%s' % (a2,a2)
            l4 = 'Missing'
            llist = [l1,l2,l3,l4]
            llist = [x.encode('ascii') for x in llist]
            plist = [pch1,pch2,pch3,pch4]
            clist = [colours[x] for x in plist]
            plist = [plotsym[x] for x in plist]
            xvec = [float(i) for i in xvec.split(',')]
            yvec = [float(i) for i in yvec.split(',')]
            rvec = [float(i) for i in rvec.split(',')]
            thetavec = [float(i) for i in thetavec.split(',')]
            zvec = [int(i) for i in pchvec.split(',')]
            colvec = [colours['%d' % i] for i in zvec]
            pchvec = [plotsym['%d' % i] for i in zvec]
            nrows = len(xvec)
            m = '%s %s chr%s:%s' % (title,marker.encode('ascii'),chrom.encode('ascii'),offset.encode('ascii'))
            ylab = 'Normalised Y'  
            xlab = 'Normalised X'
            fname = '%s_%s.pdf' % (xyfname,marker)
            pdff = os.path.join(outdir,fname)
            flist.append(fname)
            expl = '%s Normalized XY plot' % marker
            explanations.append(expl)
            rpy.r.pdf( os.path.join(outdir,fname), h , w  )
            rpy.r.par(mai=(1,1,1,0.5))
            rpy.r.plot(xvec,yvec,xlab=xlab,ylab=ylab,main=m,pch=pchvec, col=colvec)
            rpy.r.legend("topright",legend=llist,pch=plist,col=clist,title=marker)
            rpy.r.grid(nx = None, ny = None, col = "lightgray", lty = "dotted")
	    rpy.r.dev_off()
            print 'dev off'
            ylab = 'Normalised R value'  
            xlab = 'Normalised Theta value'
            fname = '%s_%s.pdf' % (rthetafname,marker)
            flist.append(fname)
            expl = '%s Theta/R plot' % marker
            explanations.append(expl)
            pdff = os.path.join(outdir,fname)
            rpy.r.pdf( pdff, h , w  )
            rpy.r.par(mai=(1,1,1,0.5))
            rpy.r.plot(thetavec,rvec,xlab=xlab,ylab=ylab,main=m,pch=pchvec, col=colvec, cex=0.8)
            rpy.r.legend("topright",legend=llist,pch=plist,col=clist,title=marker)
            rpy.r.dev_off()
        else:
            flist.append('')
            explanations.append('##ERROR: %s NOT found in full call report##' % marker)
    if mogrify:
        vcl = 'mogrify -format jpg -resize %s %s' % (mogresize, os.path.join(outdir,'*.pdf'))
        # make thumbnail images    
        x=subprocess.Popen(vcl,shell=True)
        retval = x.wait()
    return flist,explanations
    
def plotMaker():
    """
    <command interpreter="python">
        illPlotmaker.py "$rslist" $outfile1 $outfile1.files_path "$title" $height $width $mogrify      
    </command>    

        """
    progname = os.path.split(sys.argv[0])[-1] # get script name
    print >> sys.stdout,'## %s at %s cl=%s' % (progname, timenow(), ' '.join(sys.argv))
    orslist = sys.argv[1].lower().replace('__cr__',' ')
    orslist = orslist.replace('__cn__',' ').split() # galaxy replaces newlines with XX - go figure
    badmarkers = []
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
    outfile = sys.argv[2]
    newfilepath = sys.argv[3]
    title = sys.argv[4] # unchanged
    height = int(sys.argv[5])
    width = int(sys.argv[6])
    thumbnails = (sys.argv[7].lower() == "true")
    xyfname = 'X_Y'
    rthetafname = 'R_Theta'
    try:
        os.makedirs(newfilepath)
    except:
        pass
    newfiles,explanations = makeplots(mlist=rslist,title=title,outfile=outfile,outdir=newfilepath,
              xyfname=xyfname, rthetafname=rthetafname,h=height,w=width,mogrify=thumbnails)
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
    if len(sys.argv) > 1 and sys.argv[1].lower() == 'build':
        foo = illFCR() 
        foo.buildSQLite() # be warned, this takes hours!
    elif len(sys.argv) > 1 and sys.argv[1].lower() == 'buildaffy':
        foo = affyFCR(dbname='/nfs/framinghamShare_nov2007/geno/affy500raw/affy500Framingham.sqlite',
                 tablename='intensities', intensityfile='affy500Framingham.intensities',
                 gdir = '/nfs/framinghamShare_nov2007/geno/affy500raw',genoprefix='*.geno.*')
 
        foo.buildSQLite() # be warned, this takes hours!
    else:
        plotMaker()

