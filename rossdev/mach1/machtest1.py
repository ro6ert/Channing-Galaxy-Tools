#! env python
"""
 mach runner 
 many complications - training set from main to-be-inferred dataset 
 can be 300 subjects - ?founders only if family data?
 need to see what works best. What happens if trios in training data?

Assumes corresponding error/crossover file names based on incoming basename

Optionally subset the incoming ped file and fix alleles/affection for mach

If no error/crossover files:
    Optionally use whole or subset chromosome of reference haplotypes
    run mach step 1 to create error/crossover files

run mach step 2
remove low rsq markers
create ped/map file

Assume wga is like a series of 24 regions so detect and split
nsamp = 300
rounds = 50
hapsource = g1000 = '/share/shared/data/1000g/2010_06' # 2010_06 or 2010_03 or 2009_08 ?
chrom,offset,length options for subset
randfrac for subset testing

preparation:
1. If family data, write to_be_imputed as pure trios - use 'to' flag
2. If family data, filter to_be_imputed to make training ped as founders only and subset to nsamp else subset to nsamp
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
 and keep to a single chromosome or region - need to add control for a WGA input file
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
import sys,os,random,copy,tempfile,subprocess,string
from optparse import OptionParser
prog = os.path.split(sys.argv[0])[1]
debug = 1
g1000 = '/share/shared/data/1000g/2010_06' # 2010_06 or 2010_03 or 2009_08 ?


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
        mins = dict(zip(cs,[sys.maxint for x in cs]))
        maxs = dict(zip(cs,[-9 for x in cs]))
        for i,x in enumerate(clist):
            ofs = olist[i]
            if ofs < mins[x]:
                mins[x] = ofs
            if ofs > maxs[x]:
                maxs[x] = ofs  
        self.mapf.seek(0)
        self.mapl = mapl
        self.cmaxs = maxs
        self.cmins = mins # these are useful for slicing
    
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
            assert len(srow) == 6+2*self.lm,'setupPed error - short row len %d expected %d' % (len(srow),6+2*self.lm)
            ped = srow[:6]
            fid = ped[0]
            iid = ped[1]
            k = '%s_%s' % (fid,iid)
            peds.setdefault(k,copy.copy(ped))
            if ped[2] <> '0' or ped[3] <> '0': # non founder
                offspring.setdefault(k,copy.copy(ped))
                noffs.setdefault(fid,0)
                noffs[fid] += 1 # record each offspring
                if pureTrios: # optimise choice - select affected offspring if possible
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
        
    def sliceMap(self,chrom=None,offset=0,length=1000,randFrac=None):
        """ destructive slice - can be undone by calling setupMap again
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


    def writePed(self,outPath=None,subsetn=None):
        """write out a potentially filtered ped and map 
        call sliceMap first to create a sliced set of markers
        remember to reset with setupMap if needed again
        """
        assert outPath <> None,'## error subsetPed called without outpath'
        assert self.putOut in ['all','fo','oo','to'],'# error pedfile writePed called with putOut = %s - remember to call setupPed first?' % putOut
        frac = None
        if subsetn:
            assert self.putOut in ['oo','fo'],'# error in writePed - frac = %f but putOut = %s - odd?' % (frac,self.putOut)
            if self.lenped > subsetn:
                frac = float(subsetn)/self.lenped # don't subset if too few to start with
        alleledict = {'0':'0','1':'A','2':'C','3':'G','4':'T','1':'1','2':'2','3':'3','4':'4','-9':'0','N':'0','-':'0'}
        self.setupPed(putOut=putOut)
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
            if putOut == 'fo':
                if (ped[2] == '0') and (ped[3] == '0'):
                    writeme = True
            elif putOut == 'oo':
                if (ped[2] <> '0') or (ped[3] <> '0'):
                    writeme = True
            elif putOut == 'to': # trios only
                if self.trioped.get(k,None): # put out to triofile
                    writeme = True
            elif putOut == 'all':
                writeme = True
            else:
                assert False, 'Odd, putOut=%s in writePed? Nothing will be written!' % putOut
            if writeme and (frac <> None and random.random() > frac):
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
                rounds=None,inseed=0,machexe='mach1',cd=None,nmax=300,chrom=None,offset=None,length=None,**kwd):
        """
        all basenames include full paths
        basename is ped file basename
        snpsource and hapsource point to the set of reference haplotypes - hap/snps pair or if None, these will be sorted out 
        including checking to see if a subset of chrom is required
        stage1name is the name for the crossover and erate files from mach stage 1 - made if not exists
        outname is used as stage2 output name
        """
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
        assert rounds > 0,'## error machRun class instance initiated with rounds not > 0'
        assert self.basename <> None,'## error machRun class instance initiated with missing basename'
        assert os.path.exists('%s.ped' % self.basename),'## error no ped file found for machRun class instance basename %s' % basename
        assert self.refname <> None,'## error machRun class instance initiated with missing reference haplotype name'
        self.chroms,self.starts,self.ends = self.getExtent()
        assert self.chrom in chroms,'## error requested chromosome %s not in %s.map chroms (=%s) in machRun' % (self.chrom,self.basename,chroms)  
        
    def slicehaps(self,newPrefix=None):
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
        snpf = '%s/snps/chr%s.snps' % (self.g1000,self.chrom)
        assert os.path.isfile(snpf),'# error - slicehaps cannot open snp file %s' % snpf
        hapf = '%s/hap/chr%s.hap' % (self.g1000,self.chrom)
        assert os.path.isfile(hapf),'# error - slicehaps cannot open hap file %s' % hapf
        mapf = '%s/map/chr%s.map' % (self.g1000,self.chrom)
        assert os.path.isfile(mapf),'# error - slicehaps cannot open map file %s' % mapf
        mapl = open(mapf,'r').readlines() # chr name offs
        mapls = [x.split('\t') for x in mapl]
        lastone = self.offset + self.length
        outhapf = '%s.hap' % newPrefix
        outsnpf = '%s.snps' % newPrefix
        mapi = [i for i,x in enumerate(mapls) if (int(x[2]) >= self.offset and int(x[2]) <= (lastone))] # use indices
        assert len(mapi) > 0, '# error slicehaps %s has zero length from %d to %d in %s' % (mapf,offset,lastone,snpf)
        if len(mapi) == len(mapl): # not a subset - return path to whole chromosome of data
            spath = '%s/snps/chr%s.snps' % (self.g1000,self.chrom)
            hpath = '%s/hap/chr%s.hap' % (self.g1000,self.chrom)
            if os.path.exists(outhapf):
                os.system('rm -rf %s' % outhapf)
            if os.path.exists(outsnpf):
                os.system('rm -rf %s' % outsnpf)
            os.link(spath,outsnpf)
            os.link(hpath,outhapf)
        else:
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
                outhap.write(' '.join(rb)) # keep ugliness for mach
                outhap.write(' ')
                outhap.write(''.join(subh)) 
            outhap.close()
            h.close()
        self.snpsource = outsnpf
        self.hapsource = outhapf


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


    def compute(self):
        datname = '%s.dat' % self.basename
        if not os.path.exists(datname):
            self.makeDat(self.basename)
        if not os.path.exists(self.snpsource) or not os.path.exists(self.hapsource): # make
            self.slicehaps(newPrefix=self.stage1name)
        assert os.path.exists(self.snpsource),'## error no snps file found for machRun compute stage1 %s' % self.stage1name
        assert os.path.exists(self.hapsource),'## error no hap file found for machRun compute refname %s' % self.stage1name
        alog = self.doMachMain()
        return ''.join(alog)
        
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
        """remove all -9 from ped 
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

        
    def doMach1(self):
        """
        stage1 - given a few hundred independent subjects for a whole or part chromosome, corresponding regional snp and hap data 
        make crossovermap and error rates for stage 2
        
        expects a ped file at basepath

        mach1 -d ${MF}.dat -p ${MF}.ped -s ${FNAME}.snps -h ${FNAME}.hap -r $ROU -o $MF --seed $RND >& ${MF}.stage1.log
        """
        datname = '%s.dat' % self.trainname
        if not os.path.exists(datname):
            self.makeDat(self.trainname)
        seed = self.nextSeed()
        self.cleanAff(self.trainname)
        self.truncTrain(self.trainname)
        self.makeDat(self.trainname)
        cl = 'mach1 -d %s.dat -p %s.ped.short -s %s.snps -h %s.hap -r %d -o %s --seed %d --greedy' % (self.trainname,self.trainname,
          self.refname,self.refname,self.rounds,self.stage1name,seed) 
        print '## doMach1 cl = ',cl
        alog = self.run(cl)
        return alog
        
    def doMachMain(self):
        """
        need basename=None,refname=None, eratename=None, crossratename=None, outname=None,haps=None,snps=None,rounds=None,seed=0
        stage2
        mach1 -d ${MF}.dat -p ${MF}.ped -s ${FNAME}.snps -h ${FNAME}.hap --errorMap ${MF}.erate 
        --crossoverMap ${MF}.rec --greedy --autoFlip --seed $RND --mle --mldetails -o $MF >& ${MF}.impute.log 
        if no erate/crossrate files exist, make them - this could take a long time
        """
        if not (os.path.exists(self.eratename) and os.path.exists(self.crossratename)): # make
            alog = self.doMach1()
        print '## doMach1 log = ',''.join(alog)
        seed = self.nextSeed()
        self.cleanAff(self.basename)
        cl = """mach1 -d %s.dat -p %s.ped -s %s.snps -h %s.hap --errorMap %s.erate --crossoverMap %s.rec --greedy --autoFlip --seed %d --mle --mldetails -o %s""" % \
          (self.basename,self.basename,self.refname,self.refname,self.stage1name,self.stage1name,seed,self.stage2name)
        print '## doMach2 cl = ',cl
        alog = self.run(cl)
        return alog
    
class plinkRun():
    """
    plink run container - inputs, intermediates, outputs, settings 
    """   
    def __init__(self, **kwd):
        self.plinkPath = kwd.get('plinkPath','plink')
        self.baseName=kwd.get('baseName',None)
        self.baseFType=kwd.get('baseFType',None)
        self.extraParams = kwd.get('extraParams',None)
        self.outName=kwd.get('outName',None)
        self.logPath=kwd.get('outLogPath','plink.log')
        self.plinkTasks=kwd.get('plinkTasks',[])
        self.vclBase=kwd.get('vclBase','') # whatever it takes at the start of a cl to invoke plink - often 'plink'
        self.cd=kwd.get('cd','')
        self.clBase=kwd.get('clBase','')
        self.setBase = None
        self.setOut = None
        if self.baseName:
            if self.baseFType:
                if self.baseFType == 'pbed':
                    self.setBase = '--bfile %s' % self.baseName
                else:
                    self.setBase = '--file %s' % self.baseName
        if self.outName:
            self.setOut = '--out %s' % self.outName        
            
    def doThese(self,extraOptions=''):
        """generic analysis - sets up baseis    
        eg make a linkage ped file with --recode in extra options
        """
        assert self.setBase <> None,'#error plinkSettings doThese - no basePath/baseFType found'
        cl = []
        cl.append(self.setBase)
        if extraOptions > '':
            cl.append(extraOptions) # eg --recode
        if self.setOut:
            cl.append(self.setOut)
        self.plinkTasks.append(' '.join(cl))
       
    def addFilter(self,ftype='--filter-founders',extraOptions='--recode'):
        assert ftype in ['--filter-founders','--filter-nonfounders'],'#error plinkSettings addfilter - ftype = %s' % ftype
        if len(extraOptions) > 0:
            extraOptions = '%s %s' % (extraOptions,ftype)
        else:
            extraOptions = ftype
        self.doThese(extraOptions=extraOptions)
        
    def run(self):
        """
        plink blathers when doing pruning - ignore
        
        """
        if self.cd > '':
            cd = self.cd
        else:
            cd = os.getcwd()
        fplog,plog = tempfile.mkstemp()
        alog = []
        alog.append('## Rgenetics: http://rgenetics.org Galaxy Tools rgQC.py Plink runner\n')
        for task in self.plinkTasks: # each is a list
            vcl = [self.plinkPath,]
            if self.clBase > '':
                vcl.append(self.clBase)
            vcl.append(task)
            sto = file(plog,'w')
            x = subprocess.Popen(' '.join(vcl),shell=True,stdout=sto,stderr=sto,cwd=cd)
            retval = x.wait()
            sto.close()
            try:
                lplog = file(plog,'r').readlines()
                lplog = [x for x in lplog if x.find('Pruning SNP') == -1]
                alog += lplog
                alog.append('\n')
                os.unlink(plog) # no longer needed
            except:
                alog.append('### %s Strange - no std out from plink when running command line\n%s\n' % (timenow(),' '.join(vcl)))
        return alog

            

def makeFFO(basename=None):
    """create basenameFOUNDERS and basenameOFFSPRING lped files
    """
    flags = ['--filter-founders','--filter-nonfounders']
    names = ['FOUNDERS','OFFSPRING']
    extraoptions=''
    for i in range(len(names)):    
        flag = flags[i]
        outname = '%s%s' % (basename,names[i])
        p = plinkRun(baseName=basename,baseFType='lped',outName=outname)
        p.addFilter(ftype=flag,extraOptions='--recode --alleleACGT') # must recode!
        alog = p.run()
        print '#makeFFO alog1 =',''.join(alog)        
        p = plinkRun(baseName=outname,baseFType='lped',outName=outname)
        p.doThese(extraOptions='--nonfounders --freq --hardy2')
        alog = p.run()
        print '#makeFFO alog2 =',''.join(alog)

def cleanAff(basename=None):
    """remove all -9 from ped
    """
    cl = "sed 's/-9/0/g' %s.ped > %s.ped.clean" % (basename,basename)
    os.system(cl)
    os.unlink('%s.ped' % basename)
    os.rename('%s.ped.clean' % basename, '%s.ped' % basename)

def machtest():
    """eg 
    python machtest.py testCARE_trios_all test3 20000000 2000000 22 0.3 True 02445 2
    or
    /mnt/memefs/GWAS/SHARP/CARE_SHARP/CARE_SHARP

    """
    npar = 10
    if len(sys.argv) < npar:
        print 'Usage: makemachsamp.py infile_prefix outfile_prefix offset slicelength chrom randfrac onlyOne inseed rounds'
        sys.exit(1)
    prefix,outprefix,offset,length,chrom,randfrac,onlyOne,inseed,rounds = sys.argv[1:npar]
    onlyOne = onlyOne and onlyOne.lower() in ['true','1']
    if inseed == '':
        inseed = 0
    else:
        try:
            inseed = int(inseed)
        except:
            print 'inseed %s not an integer in cl = %s' % (inseed,sys.argv[1:npar])
    makeFFO(basename=outprefix) # create offspring/founders only ped files for mach phase 1
    fped = '%sFOUNDERS' % outprefix
    m = machRun(basename=outprefix,trainname=fped,refname=outprefix,stage1name=outprefix,outname=outprefix,rounds=int(rounds),seed=inseed,machexe='mach1')
    m.compute()
    
def slicetest():
    """
        python machtest.py testCARE_trios_all test3 20000000 2000000 22 0.3 True 02445 2 slicetest
    or
    /mnt/memefs/GWAS/SHARP/CARE_SHARP/CARE_SHARP

    """
    npar = 8
    if len(sys.argv) < npar:
        print 'Usage: makemachsamp.py infile_prefix outfile_prefix offset slicelength chrom randfrac onlyOne'
        sys.exit(1)
    prefix,outprefix,offset,length,chrom,randfrac,onlyOne = sys.argv[1:npar]
    onlyOne = onlyOne and onlyOne.lower() in ['true','1']
    allbase,allbasetrio = slicePed(prefix=prefix,newPrefix=outprefix,offset=int(offset),length=int(length),chrom=chrom,randfrac=None,onlyOne=True)
    allname = slicehaps(newPrefix=outprefix,chrom=chrom,offset=int(offset),length=int(length))

helptext="""
Purposes: impute WGA or a specific region using MACH
Optional: generate random subsets from a region and impute 

Takes an input ped/map, mach settings, source for reference haplotypes and optionally, 
subset chromosome start and length
-b basefile (input ped/map basename)
-t title (use as name for all new output file names - intermediate file names like crossover/erates are basename_randfrac_fo eg)
-c chrom (otherwise split into chroms and run each)
-s offset (if c)
-l length (if c)
--randfrac 0.3 (if c) default 1.0 - use eg 0.3 for testing quality
-1 fo | oo (default is use founders for stage 1 training set)
-2 all | fo | oo | puretrios default is impute all in stage 2)
-n number of fo or oo to use for stage 1 randomly chosen
-h /share/shared/data/1000g/2010-06 (source directory for haplotypes to use for imputation)
-r rounds default 50
-g use greedy option
--compact use compact option - generally not recommended according to recent docs
"""

def mrun():
    """
    """
    killme = string.punctuation + string.whitespace
    trantab = string.maketrans(killme,'_'*len(killme))
    parser = OptionParser(usage=helptext, version="%prog 0.01")
    a = parser.add_option
    a("-b","--basefile",dest="basefile", default=None,help='Ped/Map file input')
    a("-c","--chrom",dest="chrom", default=None,help='Subset to chromosome')
    a("-l","--length",dest="length",type="int",default=sys.maxint,help='Subset length')
    a("-r","--rounds",dest="rounds",type="int",default=50,help='MACH stage 1 rounds - eg 50')
    a("-s","--start",dest="offset",type="int",default=0,help='Subset start offset')
    a("-o","--outname",dest="stage2name", default='mymachjob',help='Title for output files')
    a("-p","--outpath",dest="outpath", default='.',help='Path for all output and intermediate files')    
    a("-f","--randfrac",dest="randfrac",type="float",default=1.0,help='Subset a random proportion of markers')
    a("-1","--stage1",dest="stage1",default='fo',help='Stage 1 use founders only (fo) or independent offspring only (oo) for training')
    a("-2","--stage2",dest="stage2",default='all',help='Stage 2 use all (all), founders only (fo), puretrios, or independent offspring only (oo) for imputation')
    a("-a","--happath",dest="g1000",default=None,help='Path to home of 1000 genomes or hapmap MACH compatible haplotypes to impute')
    a("-n","--nstage1",dest="nstage1",type="int",default=300,help='Stage 1 number of subjects')
    a("--greedy",dest="greedy",default=False,action='store_true',help='Add MACH greedy parameter')
    a("--compact",dest="compact",default=False,action='store_true',help='Add MACH compact parameter (not recommended)')
    a("-e","--seed",dest="seed",default=None,help='Seed for all mach operations (default is reseed and record seeds at each step)')    
    (options,args) = parser.parse_args()
    assert os.path.isfile('%s.ped' % options.basefile),'# error: %s needs a ped/map at %s' % (prog,options.basefile)
    assert os.path.isfile('%s.map' % options.basefile),'# error: %s needs a ped/map at %s' % (prog,options.basefile)
    basepath,basename = os.path.split(options.basefile) # 
    stage1outname = os.path.join(options.outpath,'%s_%s' % (basefile,'stage1'))
    stage2inname = os.path.join(options.outpath,'%s_%s_%s' % (basefile,options.stage2,'stage2in')) # will include (eg) fo
    p = pedFile(basename=options.basefile) # can now create slices
    mustsplit = False
    if len(p.cset) > 0: # multiple chromosomes
        mustsplit = True
    undoMe = False
    if chrom:
        undoMe = True
        p.sliceMap(chrom=options.chrom,offset=options.offset,length=options.length,randFrac=options.randfrac)
        stage1outname = os.path.join(options.outpath,'%s_%s_%s' % (basefile,chrom,'stage1'))
        if options.length > 0 and options.offset > 0:
            p.sliceMap(chrom=options.chrom,offset=options.offset,length=options.length,randFrac=options.randFrac)
            stage1outname = os.path.join(options.outpath,'%s_%s_%d_%d_%s' % (basefile,chrom,options.offset,options.offset+options.length,'stage1'))  
    elif mustsplit:
        print '# error - chromosome splitting not yet implemented - use the -c option to set a single chromosome please. Stopping'
        sys.exit(1)
    p.setupPed(putOut=options.stage1) 
    p.writePed(outPath=stage1outname,subsetn=options.nstage1)
    # ready for stage 1 - compute will take care of reference haplotypes
    p.setupPed(putOut=options.stage2)
    p.writePed(outPath=stage2inname)
    # ready for stage 2
    if undoMe: # remove slice and other restrictions
        p.initPed()
        p.initMap()        
    m = machRun(basename=basefile, stage1name=stage1outname, outname=title,rounds=None,inseed=int(options.seed),
        machexe='mach1',cd=None,nmax=options.nstage1,chrom=options.chrom,offset=options.offset,length=options.length,hapsource=options.happath)
    alog = m.compute()
    print alog
    
if __name__ == "__main__":    
    mrun()



