#! env python
"""
 mach runner 
 copyright July 2010 Ross Lazarus
 All rights reserved 
 For the Rgenetics project
 Released under the LGPL
 
 many complications

 python machtest.py -b tinywga -c 22 -s 21784722 -l 123000 -r 1 -o temp -p ./test -1 fo -2 fo -a /share/shared/data/1000g/2010_06 --greedy
 seems to work pretty much as expected.

 Need a training set from main to-be-inferred dataset 
 can be 300 subjects - ?founders only if family data?
 need to see what works best. What happens if trios in training data?

 Optionally subset the incoming ped file and fix alleles/affection for mach

This can create and test a random subset of markers for testing effects of density on snp imputation accuracy
ie play hide the snp.

Mach stages:

preparation:
1. If family data, can write the to_be_imputed ped as pure trios - use 'to' flag
2. If family data, filter to_be_imputed to make training ped as founders only and subset to nsamp
3. If not family data subset to_be_imputed to nsamp training set
4. If >1 chrom, split into single chrom jobs

stage 1:
if crossover and erate not available
1. Make a haps/snps set for region if is subset - otherwise use original hap/snps chrom set
2. Run mach stage 1     
otherwise recycle

stage 2:
 1. run mach stage2
 2. filter on rsq to create final ped


 has a machRun class to hide some details
 and keep to a single chromosome or region - control for a WGA input file not yet tested but done
 includes code to obtain cut down hap/snps files if region and a
 generic way to cut down a ped file
 including optionally taking a random sample of markers over an inteval
 also includes a new super clumsy plink runner - still not sure it isn't easier to just bang command lines together 
 and run them...
 started out life as a
 script to make some small MACH samples
 both family and independent to test differences in quality
 need slice of /udd/relli/1000G_Sanger_0908/0908_CEU_NoSingleton/hap/0908_CEU_NoSingleton_chr22.hap eg
 actually, updated to q:/share/shared/data/1000g
 to match a full slice from /mnt/memefs/GWAS/CAMP_610/CAMP_610.ped
 then a random subset from the camp slice to run mach and compare the imputed with genotyped markers

 use with runme.sh to generate some frq, hwe etc from subsets - founders, non-founders, all
"""
import sys,os,random,copy,tempfile,subprocess,string,time,logging,datetime,gzip
from optparse import OptionParser
prog = os.path.split(sys.argv[0])[1]
debug = 1
g1000 = '/share/shared/data/1000g/2010_06' # 2010_06 or 2010_03 or 2009_08 ?

def setLogging(appname=prog,outdir='.'):
    """setup a logger
    """
    logdir = outdir
    if not os.path.isdir(logdir):
        os.makedirs(logdir)
    today=datetime.datetime.now().strftime('%Y%m%d')
    logging.basicConfig(level=logging.INFO,
                format='%(asctime)s %(levelname)s %(message)s',
                datefmt='%a, %d %b %Y %H:%M:%S',
                filename=os.path.join(logdir,'%s_%s.log' % (appname,today)),
                filemode='a')
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    formatter = logging.Formatter('%(levelname)-4s %(message)s')
    console.setFormatter(formatter)
    logging.getLogger('').addHandler(console)


def timenow():
    """return current time as a string
    """
    return time.strftime('%d/%m/%Y %H:%M:%S', time.localtime(time.time()))

def fTranspose(inf='',ofname='',insep=' ',outsep=' '):
    """quick and dirty transpose - not needed if we only use mlgeno but for mlinfo dose munging
    all in ram - can't use column = file trickery
    """
    assert os.path.isfile(inf),'# Error - fTranspose cannot open supplied infile = %s' % inf
    f = open(inf,'r')
    d = []
    rowlen = None
    logging.info('fTranspose reading %s' % inf)
    for i,s in enumerate(f):
        row = s.strip().split(insep)
        if i == 0: # first row
           rowlen = len(row)
        else:
            assert rowlen == len(row),'# Error - fTranspose inf %s mixed rows %d; at %d, %d' % (inf,rowlen,i,len(row)) 
        d.append(row)
    logging.info('transposing')
    d = zip(*d)
    logging.info( 'writing %s' % ofname)
    fo = open(ofname,'w')
    for row in d:
      fo.write(outsep.join(row))
      fo.write('\n')
    fo.close()

class pedFile():
    """ here we go again
    specialised slicer - plink cannot currently autoprune a file to give optimal complete trios
    eg to slice a random 50% fraction from chr22 only and write founders only
    p = pedFile(basename='/opt/galaxy/test-data/tinywga')
    p.subsetPed(chrom='22',offset=0,length=999999999,randFrac=0.5,putOut='fo',outPath='/tmp/foo')    
    """
    def __init__(self,basename=None,**kwd):
        """
        """
        self.basename = basename
        self.mapname = '%s.map' % basename
        self.pedname = '%s.ped' % basename
        assert os.path.isfile(self.mapname),'## error pedfile instantiated but no map available at %s' % mapfile
        assert os.path.isfile(self.pedname),'## error pedfile instantiated but no ped available at %s' % pedfile
        self.pfile = open(self.pedname,'r')
        self.mapf = open(self.mapname,'r') 
        self.initMap()
        self.initPed()
        
    def initMap(self):
        """prepare for slice"""
        mapl = self.mapf.readlines()
        self.lmap = [x.strip().split() for x in mapl if len(x.split()) >=4]
        self.rslist = [x[1] for x in self.lmap]
        self.clist = [x[0] for x in self.lmap]
        self.cset = set(self.clist)
        self.olist = [x[3] for x in self.lmap]
        self.mslice = range(len(self.rslist)) # index for each active marker - can be changed by slice
        self.lenmap = len(self.lmap)
        mins = dict(zip(self.cset,[sys.maxint for x in self.cset]))
        maxs = dict(zip(self.cset,[-9 for x in self.cset]))
        for i,x in enumerate(self.clist):
            ofs = self.olist[i]
            if ofs < mins[x]:
                mins[x] = ofs
            if ofs > maxs[x]:
                maxs[x] = ofs  
        self.mapf.seek(0)
        self.mapl = mapl
        self.cmaxs = maxs
        self.cmins = mins # these are useful for slicing
        
    def sliceMap(self,chrom=None,offset=0,length=sys.maxint,randFrac=None):
        """ destructive slice - can be undone by calling sliceMap again
        """
        assert chrom in self.cset,'# error sliceMap called with chrom %s not found in cset %s' % (chrom,self.cset)
        end = offset + length
        mslice = [i for i,x in enumerate(self.lmap) if (x[0]==chrom and int(x[3]) >= offset and int(x[3]) <= (end))] 
        if randFrac > 0 and randFrac < 1.0: # random subset of all markers for testing eg
            n = int(len(mslice)*randFrac)
            assert n > 0,'sliceMap error %f of %d < 1' % (randFrac,len(mslice))
            randm = random.sample(mslice,n) # subset of marker indices
            randm.sort() # sample should retain order but..doesn't seem to?
            mslice = randm
        assert len(self.mslice) > 0,'# error sliceMap created an empty new map!'
        self.mslice = mslice
        self.lenmap = len(mslice)
        self.rslist = [self.rslist[i] for i in mslice]
        self.clist = [self.clist[i] for i in mslice]
        self.olist = [self.olist[i] for i in mslice] # rebuild indices
    
    def initPed(self):
        """ prepare pedigree structures needed for setupPed which
        may need to be run multiple times per pedFile
        """       
        peds = {}
        offspring = {}
        trioped = {}
        noffs = {}
        affofs = {}
        for row in self.pfile: # need > 1 pass to find best (one of the affected offspring) choices for strict trios
            srow = row.split()
            assert len(srow) == 6+2*self.lenmap,'setupPed error - short row len %d expected %d' % (len(srow),6+2*self.lenmap)
            ped = srow[:6]
            fid = ped[0]
            iid = ped[1]
            k = '%s_%s' % (fid,iid)
            peds.setdefault(k,copy.copy(ped))
            if ped[2] <> '0' or ped[3] <> '0': # non founder
                offspring.setdefault(k,copy.copy(ped))
                noffs.setdefault(fid,0)
                noffs[fid] += 1 # record each offspring
               # optimise choice - select affected offspring if possible
                if ped[5] == '2': # affected
                    affofs.setdefault(fid,[])
                    affofs[fid].append(iid)
        self.pfile.seek(0)
        self.hastrios = False
        if len(offspring) > 0:
            self.hastrios = True
        self.peds = peds
        self.offspring = offspring
        self.trioped = trioped
        self.noffs = noffs
        self.affofs = affofs
        self.lenped = len(peds)

    def setupPed(self,putOut='all'):
        """prepare for subset"""
        pureTrios = False
        if putOut == 'to':
            pureTrios = True
        fidDone = {}  
        self.putOut = putOut          
        for ped in self.offspring.values(): # flag offspring for inclusion or not
            fid,iid,pid,mid = ped[:4]
            dadk = '%s_%s' % (fid,pid)
            mumk = '%s_%s' % (fid,mid)
            kidk = '%s_%s' % (fid,iid)
            dad = self.peds.get(dadk,None)
            mum = self.peds.get(mumk,None)
            kid = self.peds.get(kidk,None)        
            if ((dad != None) or (mum != None)) and (kid != None):
                saveme = 1 # save all offspring if non onlyOne
                if pureTrios:
                    if ((dad == None) or (mum == None)):
                        saveme = 0 # don't save non nuclear trios
                    if fidDone.get(fid,None):
                        saveme = 0 # only one per family
                    else:
                        a = self.affofs.get(fid,None) # see if there are any
                        if a <> None:
                            if iid not in a: # there are affected and this is not one don't write
                                saveme = 0
                if saveme:
                    fidDone.setdefault(fid,1) # flag no more
                    self.trioped.setdefault(dadk,1) # include this complete trio in _trio file
                    self.trioped.setdefault(mumk,1)
                    self.trioped.setdefault(kidk,1)
        

    def writePed(self,outPath=None,subsetn=None):
        """write out a potentially filtered ped and map 
        remember to call sliceMap to set up the markers to be written out
        and setupPed to set up the pedigree - eg for fo or puretrios to be imputed
        """
        assert outPath <> None,'## error writePed called without outpath'
        assert self.putOut in ['all','fo','oo','to'],\
          '# error pedfile writePed called with putOut = %s - remember to call setupPed first?' % putOut
        frac = None
        if subsetn:
            assert self.putOut in ['oo','fo'],'# error in writePed - frac = %f but putOut = %s - odd?' % (frac,self.putOut)
            if self.lenped > subsetn:
                frac = float(subsetn)/self.lenped # don't subset if too few to start with
        alleledict = {'0':'0','1':'A','2':'C','3':'G','4':'T','A':'A','C':'C','G':'G','T':'T','-9':'0','N':'0','-':'0'}
        outpedf = open('%s.ped' % outPath,'w')
        outmapf = open('%s.map' % outPath,'w')
        newmap = [self.mapl[x] for x in self.mslice] # haha
        outmapf.write(''.join(newmap))
        outmapf.write('\n')
        outmapf.close()
        cols = [6+2*x for x in self.mslice]
        cols += [7+2*x for x in self.mslice]
        cols.sort()
        self.pfile.seek(0)
        recode = None
        for row in self.pfile: # finally can write out what we are looking for
            srow = row.split()
            ped = srow[:6]
            k = '%s_%s' % (ped[0],ped[1])
            writeme = False
            if self.putOut == 'fo':
                if (ped[2] == '0') and (ped[3] == '0'):
                    writeme = True
            elif self.putOut == 'oo':
                if (ped[2] <> '0') or (ped[3] <> '0'):
                    writeme = True
            elif self.putOut == 'to': # trios onlyself
                if self.trioped.get(k,None): # put out to triofile
                    writeme = True
            elif self.putOut == 'all':
                writeme = True
            else:
                assert False, 'Odd, putOut=%s in writePed? Nothing will be written!' % self.putOut
            if writeme and frac <> None: # random sample to density
                writeme = (random.random() > frac)
            if writeme:
                g = [srow[x] for x in cols]
                if recode == None:
                    try:
                        test = [int(x) for x in g[:100]] # are there number alleles?
                    except:
                        test = []
                    if len(test) > 0:
                        recode = True
                    else:
                        recode = False
                if recode:
                    g = [alleledict.get(x.upper(),'?') for x in g]  
                if ped[5] == '-9':
                    ped[5] = '0' # mach does not like -9 for affection
                ped += g                
                outpedf.write(' '.join(ped))
                outpedf.write('\n')
        outpedf.close()

class machRun():
    """
    mach run container - inputs, intermediates, outputs, settings for a chromosome or part thereof
    """  
    def __init__(self, basename=None, stage1name=None, snpsource=None, hapsource=None, outname=None,
                stage2inname=None,rounds=None,inseed=0,machexe='mach1',cd=None,nmax=300,mask=None,
                chrom=None,offset=None,length=None,outpath=None,**kwd):
        """
        all basenames include full paths
        basename is ped file basename
        snpsource and hapsource point to the set of reference haplotypes - hap/snps pair or if None, these will be sorted out 
        including checking to see if a subset of chrom is required
        stage1name is the name for the crossover and erate files from mach stage 1 - made if not exists
        outname is used as stage2 output name
        """
        self.outpath = outpath
        self.mask = mask
        self.stage2inname = stage2inname
        self.hapsource=hapsource
        self.chrom = chrom
        self.offset = offset
        self.length = length
        self.cd = cd
        self.nmax = nmax
        self.machexe = machexe
        self.basename = basename     
        self.hapsource = hapsource # reference haps for imputation
        self.snpsource = snpsource # reference snp names for imputation
        self.rounds = rounds
        self.inseed = inseed
        self.stage1name = stage1name
        self.outname = outname
        self.crossratename = stage1name
        self.eratename = stage1name
        self.greedy = kwd.get('greedy',False)
        self.compact = kwd.get('compact',False)
        assert rounds > 0,'## error machRun class instance initiated with rounds not > 0'
        assert self.basename <> None,'## error machRun class instance initiated with missing basename'
        assert os.path.exists('%s.ped' % self.basename),'## error no ped file found for machRun class instance basename %s' % basename
        self.chroms,self.starts,self.ends = self.getExtent()
        assert self.chrom in self.chroms,'## error requested chromosome %s not in %s.map chroms (=%s) in machRun' % (self.chrom,self.basename,self.chroms)  

        
    def slicehaps(self,newPrefix=None,mapPrefix=None):
        """ specific to mach hap format - take a slice from 1000g data corresponding to our slice
        or link to source haps if the slice is the entire chromosome
        note up to 2 fields may preceed the haplotype string !
        Sun,Jul 18 at 12:13am head /share/shared/data/1000g/map/chr22.map
        22	chr22:14431249	14431249
        22	rs6518413	14432239
        22	chr22:14433730	14433730
        Sun,Jul 18 at 12:13am head /share/shared/data/1000g/snps/chr22.snps
        chr22:14431249
        rs6518413
        chr22:14433730
        rs2844885
        chr22:14435122
        """
        assert newPrefix <> None, '# error - slicehaps needs a prefix for output'
        snpf = '%s/snps/chr%s.snps' % (self.hapsource,self.chrom)
        assert os.path.isfile(snpf),'# error - slicehaps cannot open snp file %s' % snpf
        hapf = '%s/hap/chr%s.hap' % (self.hapsource,self.chrom)
        assert os.path.isfile(hapf),'# error - slicehaps cannot open hap file %s' % hapf
        mapf = '%s/map/chr%s.map' % (self.hapsource,self.chrom)
        assert os.path.isfile(mapf),'# error - slicehaps cannot open map file %s' % mapf
        mapl = open(mapf,'r').readlines() # chr name offs
        mapls = [x.split('\t') for x in mapl]
        lastone = self.offset + self.length
        outhapf = '%s_ref.hap' % newPrefix
        outsnpf = '%s_ref.snps' % newPrefix
        outmapf = '%s_ref.map' % mapPrefix
        mapi = [i for i,x in enumerate(mapls) if (int(x[2]) >= self.offset and int(x[2]) <= (lastone))] # use indices
        assert len(mapi) > 0, '# error slicehaps %s has zero length from %d to %d in %s' % (mapf,offset,lastone,snpf)
        if len(mapi) == len(mapl): # not a subset - return path to whole chromosome of data
            if os.path.exists(outhapf):
                os.system('rm -rf %s' % outhapf)
            if os.path.exists(outmapf):
                os.system('rm -rf %s' % outmapf)
            if os.path.exists(outsnpf):
                os.system('rm -rf %s' % outsnpf)
            os.link(snpf,outsnpf)
            os.link(hapf,outhapf)
            os.link(mapf,outmapf)
        else:
            maps = open(mapf,'r').readlines()
            outmaps = [maps[i] for i in mapi]
            outmap = open(outmapf,'w')
            outmap.write(''.join(outmaps)) # create new mach haplo map
            snps = open(snpf,'r').readlines()
            outsnps = [snps[i] for i in mapi] # subset	x
            outsnp = open(outsnpf,'w')
            outsnp.write(''.join(outsnps))
            outsnp.close()
            outhap = open(outhapf,'w')
            h = open(hapf,'r')
            for row in h: # this is not as easy as you might think
                r = row.strip().split() # may be up to 3 - take last
                rh = r[-1] # hap
                rb = r[:-1] # start
                subh = [rh[i] for i in mapi] # ugh - but that's the format - optional row id(s)
                subh.append('\n')
                outhap.write('\t'.join(rb)) # keep ugliness for mach - note MUST be tab
                outhap.write('\t')
                outhap.write(''.join(subh)) 
            outhap.close()
            h.close()
        return outsnpf,outhapf,outmapf
    
    def doMach1(self):
        """
        stage1 - given a few hundred independent subjects for a whole or part chromosome, corresponding regional snp and hap data 
        make crossovermap and error rates for stage 2
    
        """
        datname = '%s.dat' % self.stage1name
        if not os.path.exists(datname):
            self.makeDat(self.stage1name)
        seed = self.nextSeed()
        cl = 'mach1 -d %s.dat -p %s.ped -s %s -h %s -o %s --seed %d -r %d' % (self.stage1name,self.stage1name,
            self.snpsource,self.hapsource,self.stage1name,seed,self.rounds) 
        if self.compact:
            cl = '%s --compact' % cl
        if self.greedy:
            cl = '%s --greedy' % cl
        logging.info( '## doMach1 cl = %s' % cl)
        alog = self.run(cl)
        return alog
        
    def compute(self):
        """
        need basename=None,refname=None, eratename=None, crossratename=None, outname=None,haps=None,snps=None,rounds=None,seed=0
        stage2
        mach1 -d ${MF}.dat -p ${MF}.ped -s ${FNAME}.snps -h ${FNAME}.hap --errorMap ${MF}.erate 
        --crossoverMap ${MF}.rec --greedy --autoFlip --seed $RND --mle --mldetails -o $MF >& ${MF}.impute.log 
        if no erate/crossrate files exist, make them - this could take a long time
        """
        datname = '%s.dat' % self.basename
        if not os.path.exists(datname):
            self.makeDat(self.basename)
        if self.snpsource == None or (not os.path.exists(self.snpsource) or not os.path.exists(self.hapsource)):
            outsnpf,outhapf,outmapf = self.slicehaps(newPrefix=self.stage1name,mapPrefix=self.outname) 
            # reference haplotypes - eg 1000g project for mach
            self.snpsource = outsnpf
            self.hapsource = outhapf
            self.mapsource = outmapf
        assert os.path.exists(self.snpsource),'## error no ref snps file found for machRun compute stage1 %s' % self.stage1name
        assert os.path.exists(self.hapsource),'## error no ref hap file found for machRun compute refname %s' % self.stage1name
        assert os.path.exists(self.mapsource),'## error no ref map file found for machRun compute refname %s' % self.stage1name
        alog1 = []
        if not (os.path.exists(self.eratename) and os.path.exists(self.crossratename)): # make
            alog1 = self.doMach1()
        datname = '%s.dat' % self.stage2inname
        if not os.path.exists(datname):
            self.makeDat(self.stage2inname)
        seed = self.nextSeed()
        cl = """mach1 -d %s.dat -p %s.ped -h %s -s %s --errorMap %s.erate --crossoverMap %s.rec --autoFlip --seed %d --mle --mldetails -o %s""" % \
          (self.stage2inname,self.stage2inname,self.hapsource,self.snpsource,self.stage1name,self.stage1name,seed,self.outname)
        if self.compact:
            cl = '%s --compact' % cl
        if self.greedy:
            cl = '%s --greedy' % cl
        if self.mask <> None:
            cl = '%s --mask %f' % (cl,self.mask)
        logging.info( '## doMach2 cl = %s' % cl)
        alog = self.run(cl)
        if alog1 > []:
            alog1 += alog
            return alog1
        else:
            return alog

        
    def run(self,cl=''):
        """        
        """
        if self.cd > '':
            cd = self.cd
        else:
            cd = os.getcwd()
        fplog,plog = tempfile.mkstemp()
        alog = []
        alog.append('## Rgenetics: http://rgenetics.org Galaxy Tools MACH runner\n')
        sto = file(plog,'w')
        x = subprocess.Popen(cl,shell=True,stdout=sto,stderr=sto,cwd=cd)
        retval = x.wait()
        sto.close()
        try:
            lplog = file(plog,'r').readlines()
            alog += lplog
            alog.append('\n')
            os.unlink(plog) # no longer needed
        except:
            alog.append('### %s Strange - no std out from MACH when running command line\n%s\n' % (timenow(),' '.join(vcl)))
        return alog
        
    def nextSeed(self):
        if self.inseed <> 0:
            return self.inseed
        else:
            return random.randint(1,32767)
            
    def makeDat(self,fname=None):
        """ make a dumb old .dat file for mach to match the map - note 
        additional affection pheno field so mach can read plink style 6 column pedigrees
        """
        assert fname <> None, '# error makeDat called with None as fname'
        mapname = '%s.map' % fname
        assert os.path.isfile(mapname) or os.path.islink(mapname)
        m = open(mapname,'r').readlines()
        mlist = [x.strip().split() for x in m if len(x.strip().split()) >=3]
        rslist = ['M %s' % x[1] for x in mlist]
        datname = '%s.dat' % fname
        os.system('rm -rf %s' % datname)
        outf = file(datname,'w')
        outf.write('A Affection\n')
        outf.write('\n'.join(rslist))
        outf.write('\n')
        outf.close()        
        
    def cleanAff(self,fname=None):
        """remove all -9 from ped - this is now done in the pedFile class
        """
        assert fname <> None, '# error cleanAff called with None as fname'
        cl = "sed 's/-9/0/g' %s.ped > %s.ped.clean" % (fname,fname)
        os.system(cl)
        os.unlink('%s.ped' % fname)
        os.rename('%s.ped.clean' % fname,'%s.ped' % fname)
        
    def truncTrain(self,fname=None):
        """and truncate to (say) 200 rows
        """
        cl = 'head -%d %s.ped > %s.ped.short' % (self.nmax,fname,fname)
        os.system(cl)

    def getExtent(self):
        """ find bounds of mapfile"""
        mapname = '%s.map' % self.basename
        maps = open(mapname,'r').readlines()
        m = [x.strip().split() for x in maps if len(x.strip().split()) >= 4]
        o = [int(x[3]) for x in m] # all offsets
        c = [x[0] for x in m] # all chroms
        cs = set(c)
        mins = dict(zip(cs,[sys.maxint for x in cs]))
        maxs = dict(zip(cs,[0 for x in cs]))
        for i,x in enumerate(c):
            ofs = o[i]
            if ofs < mins[x]:
                mins[x] = ofs
            if ofs > maxs[x]:
                maxs[x] = ofs          
        return cs, [mins[x] for x in cs], [maxs[x] for x in cs]

        
        
def mlinfoConv(inmachfilepath=None,outpath=None,r2thresh=0.3,qualthresh=0.0,dped={}):
    """ 
    simplest approach - use mlinfo to filter reasonable rsq values for the 'best guess' alleles
    probably easier to do it ourselves than futz with gengen perl crap
    Thu,Jul 22 at 12:04am head /share/shared/data/1000g/2010_06/map/chr22.map 
    22      chr22:14431249  14431249
    22      rs6518413       14432239
    22      chr22:14433730  14433730
    Sat,May 22 at 12:12pm head *chr22.mlinfo
    SNP     Al1     Al2     Freq1   MAF     Quality Rsq
    rs11089130      C       G       0.2653  0.2653  0.5410  0.0176
    rs738829        A       G       0.4705  0.4705  0.5037  0.0549
    rs915674        A       G       0.0736  0.0736  0.8583  0.0050
    """
    assert inmachfilepath <> None,'# error mlinfoConv given None as inmachfilepath'
    logging.info('## mlinfoConv doing %s' % inmachfilepath)
    base_name = os.path.basename(inmachfilepath)
    inmap = open('%s_ref.map' % inmachfilepath,'r').readlines()
    mli = open('%s.mlinfo' % inmachfilepath,'r').readlines()
    mli = [x.split() for x in mli if len(x.split()) > 0]
    logging.info('len dped = %d; len map = %d; len mlinfo = %d' % (len(dped),len(inmap),len(mli)))
    mlih = mli[0]
    mlih = [x.lower() for x in mlih]
    rsqpos = mlih.index('rsq')
    qualpos = mlih.index('quality')
    mli = mli[1:] # drop header
    passing = []
    failing = []
    useme = [] # indices
    for j,row in enumerate(mli):
        if float(row[rsqpos]) >= r2thresh and float(row[qualpos]) >= qualthresh:
           passing.append(row)
           useme.append(j) # index
        else:
           failing.append(row)
    outroot = os.path.join(outpath,base_name)
    outpedfroot = '%s_clean' % outroot
    outped = open('%s.ped' % outpedfroot,'w')
    outmap = open('%s.map' % outpedfroot,'w')
    logging.info('## len(useme)=%d len(inmap)=%d' % (len(useme),len(inmap)))
    for j in useme:
        row = inmap[j].split() # use this one
        row.insert(2,'0') # fake a genetic map 0 to make this plink map style
        outmap.write('\t'.join(row))
        outmap.write('\n')
    outmap.close()
    ingeno = '%s.mlgeno' % inmachfilepath
    ingenogz = '%s.mlgeno.gz' % inmachfilepath
    if os.path.exists(ingenogz):
        mlg = gzip.open(ingenogz)
    else:
        mlg = open(ingeno,'r') # big
    for pedi,row in enumerate(mlg): # 500003->20100003 ML_GENO G/G A/G G/G C/C
        row = row.split()
        if len(row) > 0:
           genos = []
           iid,fid = row[0].split('->')
           k = '%s_%s' % (iid,fid)
           thisped = dped.get(k,None) # 
           assert thisped <> None, '# error - ID %s not found in dped' % k
           assert row[1] == 'ML_GENO','## file %s row %d col 2 %s <> ML_GENO' % ('%s.mlgeno' % base_name,i,row[1])
           row = row[2:]
           for j in useme:
                g = row[j] # keep indices in useme
                g1,g2 = g.split('/')
                genos.append(g1)
                genos.append(g2)
           s = ' '.join(thisped)
           outped.write('%s %s\n' % (s,' '.join(genos)))
    outped.close()
    return outpedfroot

def machToPlink(chroms={},inbase=None,plinke='plink',outpath=None,r2thresh=0.3,qualthresh=0.0,dped={}):
    """
    some code to finalize mach imputation outputs
    """
    fp,basename = os.path.split(inbase) # get basename
    outpedroots = []
    for chrom in chroms:
        inname = os.path.join(outpath,'%s_chr%s' % (basename,chrom))
        outpedfroot = mlinfoConv(inmachfilepath=inname,outpath=outpath,r2thresh=r2thresh,
                 qualthresh=qualthresh,dped=dped)
        outpedroots.append(os.path.abspath(outpedfroot))

    # now amalgamate
    #x,allfiles = tempfile.mkstemp(prefix='plinkallfiles')
    paf = 'plinkmerge_%s.txt' % basename
    plinklist = ['%s.ped %s.map' % (x,x) for x in outpedroots[1:]] 
    # leave the first one out - this is what the merge-list option wants
    f = open(os.path.join(outpath,paf),'w')
    f.write('\n'.join(plinklist))
    f.write('\n')
    f.close()
    cl = '%s --noweb --file %s --make-bed --merge-list %s --out %s_clean' % (plinke,outpedroots[0],paf,inbase)
    logging.info( '## machtoplink calling %s' % cl)
    p = subprocess.Popen(cl,shell=True,cwd=outpath)
    retval = p.wait() # run plink
    logging.info('## machtoplink done')
        
helptext="""This is a wrapper for MACH - see http://genome.sph.umich.edu/wiki/MaCH
Purposes: impute WGA or a specific region using MACH
Optional: generate random subsets from a region and impute 

Takes an input ped/map, mach settings, source for reference haplotypes and optionally, 
subset chromosome start and length
compact option - generally not recommended according to recent docs
"""

def mrun():
    """mach runner    
    """
    def runOneChrom(chrom=None,p=None,mask=None):
        """ broken out to make iterating over multiple chromosomes simpler
        """
        stage1outname = os.path.join(options.outpath,'%s_chr%s_%s' % (options.outname,chrom,'stage1_in'))
        stage2outname = os.path.join(options.outpath,'%s_chr%s' % (options.outname,chrom))
        p.setupPed(putOut=options.stage1) 
        p.writePed(outPath=stage1outname,subsetn=options.nstage1)
        # ready for stage 1 - compute will take care of reference haplotypes
        p.setupPed(putOut=options.stage2)
        p.writePed(outPath=stage2outname)
        # ready for stage 2
        m = machRun(basename=options.basefile, stage1name=stage1outname, outname=stage2outname,
            rounds=options.rounds,inseed=int(options.seed),mask=mask,
            machexe=machexe,cd=None,nmax=options.nstage1,chrom=chrom,offset=options.offset,
            length=options.length,hapsource=options.happath,
            stage2inname=stage2outname,outpath=options.outpath)
        alog = m.compute()
        logging.info(''.join(alog))


    killme = string.punctuation + string.whitespace
    trantab = string.maketrans(killme,'_'*len(killme))
    parser = OptionParser(usage=helptext, version="%prog 0.01")
    if len(sys.argv) == 1:
        sys.argv.append('-h')
    a = parser.add_option
    a("--machexe",dest="machexe", default='mach1',help='Path to MACH executable')
    a("-b","--basefile",dest="basefile", default=None,help='Ped/Map file input')
    a("-m","--mask",dest="mask", default=None,type='float',help='MACH mask parameter for testing accuracy')
    a("-c","--chrom",dest="chrom", default=None,help='Subset to chromosome')
    a("-l","--length",dest="length",type="int",default=-1,help='Subset length')
    a("-r","--rounds",dest="rounds",type="int",default=50,help='MACH stage 1 rounds - eg 50')
    a("-s","--start",dest="offset",type="int",default=-1,help='Subset start offset')
    a("-o","--outname",dest="outname", default='mymachjob',help='Title for output files')
    a("-p","--outpath",dest="outpath", default='.',help='Path for all output and intermediate files')    
    a("-f","--randfrac",dest="randfrac",type="float",default=1.0,help='Subset a random proportion of markers')
    a("-1","--stage1",dest="stage1",default='fo',
      help='Stage 1 use founders only (fo) or independent offspring only (oo) for training')
    a("-2","--stage2",dest="stage2",default='all',
      help='Stage 2 use all (all), founders (fo), puretrios, or offspring only (oo) for imputation')
    a("-a","--happath",dest="happath",default=None,
      help='Path to home of 1000 genomes or hapmap MACH compatible haplotypes to impute')
    a("-n","--nstage1",dest="nstage1",type="int",default=300,help='Stage 1 number of subjects')
    a("--greedy",dest="greedy",default=False,action='store_true',help='Add MACH greedy parameter')
    a("--compact",dest="compact",default=False,action='store_true',
      help='Add MACH compact parameter (not recommended)')
    a("-e","--seed",dest="seed",default=-1,type="int",
      help='Seed for all mach operations (default is reseed and record seeds at each step)')    
    (options,args) = parser.parse_args()
    assert options.basefile <> None,'# error: %s must have an input ped/map basename supplied %s' % (prog,options.basefile)
    assert os.path.isfile('%s.ped' % options.basefile),'# error: %s needs a ped/map at %s' % (prog,options.basefile)
    assert os.path.isfile('%s.map' % options.basefile),'# error: %s needs a ped/map at %s' % (prog,options.basefile)
    if not os.path.exists(options.outpath):
        os.makedirs(options.outpath)
    setLogging(appname=options.outname,outdir=options.outpath)
    basepath,basename = os.path.split(options.basefile) # 
    p = pedFile(basename=options.basefile) # can now create slices
    if len(p.cset) == 1 and options.chrom == None: # input has only one
        options.chrom = p.cset[0]
    if options.chrom <> None:
        p.sliceMap(chrom=options.chrom,offset=options.offset,length=options.length,randFrac=options.randfrac)
        runOneChrom(chrom=options.chrom,p=p,mask=options.mask) 
    else:
        logging.warn('# chromosome splitting - this will take a while')
        for c in p.cset: # for each chrom we found in the pedfile - assume all are full chromosomes?
            p.sliceMap(chrom=options.chrom,randFrac=options.randfrac) # use default offset/length = all
            runOneChrom(chrom=c,p=p,mask=options.mask)    
    machToPlink(chroms=p.cset,inbase=options.outname,plinke='plink',outpath=options.outpath,r2thresh=0.3,qualthresh=0.0,dped=p.peds)
    logging.shutdown()
    
if __name__ == "__main__":    
    """
    python machtest.py -b tinywga -c 22 -s 21784722 -l 123000 -r 1 -o temp -p ./test -1 fo -2 fo -a /share/shared/data/1000g/2010_06 --greedy
    produces:
    Thu,Jul 22 at 1:29am ls -l test
total 227
-rw-r--r-- 1 rerla re 11524 Jul 22 00:54 alltest_20100722.log
-rw-r--r-- 1 rerla re  5171 Jul 22 00:54 alltest_chr22_clean.map
-rw-r--r-- 1 rerla re 34363 Jul 22 00:54 alltest_chr22_clean.ped
-rw-r--r-- 1 rerla re   312 Jul 22 00:35 alltest_chr22.dat
-rw-r--r-- 1 rerla re  6288 Jul 22 00:54 alltest_chr22.erate
-rw-r--r-- 1 rerla re   601 Jul 22 00:54 alltest_chr22.map
-rw-r--r-- 1 rerla re 17908 Jul 22 00:54 alltest_chr22.mldose.gz
-rw-r--r-- 1 rerla re  3222 Jul 22 00:54 alltest_chr22.mlgeno.gz
-rw-r--r-- 1 rerla re 10798 Jul 22 00:54 alltest_chr22.mlinfo
-rw-r--r-- 1 rerla re 29512 Jul 22 00:54 alltest_chr22.mlprob.gz
-rw-r--r-- 1 rerla re 15636 Jul 22 00:54 alltest_chr22.mlqc.gz
-rw-r--r-- 1 rerla re  4603 Jul 22 00:54 alltest_chr22.ped
-rw-r--r-- 1 rerla re  5290 Jul 22 00:54 alltest_chr22.rec
-rw-r--r-- 1 rerla re  5764 Jul 22 00:54 alltest_chr22_ref.map
-rw-r--r-- 1 rerla re   312 Jul 22 00:35 alltest_chr22_stage1_in.dat
-rw-r--r-- 1 rerla re  6288 Jul 22 00:54 alltest_chr22_stage1_in.erate
-rw-r--r-- 1 rerla re   601 Jul 22 00:54 alltest_chr22_stage1_in.map
-rw-r--r-- 1 rerla re  3103 Jul 22 00:54 alltest_chr22_stage1_in.ped
-rw-r--r-- 1 rerla re  5290 Jul 22 00:54 alltest_chr22_stage1_in.rec
-rw-r--r-- 1 rerla re 31080 Jul 22 00:54 alltest_chr22_stage1_in_ref.hap
-rw-r--r-- 1 rerla re  2764 Jul 22 00:54 alltest_chr22_stage1_in_ref.snps
-rw-r--r-- 1 rerla re  2113 Jul 22 00:54 alltest_clean.bed
-rw-r--r-- 1 rerla re  6015 Jul 22 00:54 alltest_clean.bim
-rw-r--r-- 1 rerla re   603 Jul 22 00:54 alltest_clean.fam
-rw-r--r-- 1 rerla re  1973 Jul 22 00:54 alltest_clean.log
-rw-r--r-- 1 rerla re     1 Jul 22 00:54 plinkmerge_alltest.txt

    """
    mrun()



