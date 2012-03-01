#! env python
# script to make some small MACH samples
# both family and independent to test differences in quality
# need slice of /udd/relli/1000G_Sanger_0908/0908_CEU_NoSingleton/hap/0908_CEU_NoSingleton_chr22.hap eg
# actually, updated to q:/share/shared/data/1000g
# to match a full slice from /mnt/memefs/GWAS/CAMP_610/CAMP_610.ped
# then a random subset from the camp slice to run mach and compare the imputed with genotyped markers

# use with runme.sh to generate some frq, hwe etc from subsets - founders, non-founders, all
"""
#!/bin/bash
# RANDOM requires bash
DEBUG=1
#experimental subset for mach family vs independent testing
#this selects 358 markers from chr22 care_sharp and then makes a random subset of 40% of the markers
ROU=5
#for quick testing
DOALL=0
OFILE=test2
CARE=/mnt/memefs/GWAS/SHARP/CARE_SHARP/CARE_SHARP
#otherwise use
#was python makemachsamp.py $CARE testCARE 20000000 2000000 22 0.3
#python makemachsamp.py testCARE_all $OFILE 20000000 2000000 22 0.3
#this resamples the testCARE_all at 0.3 density for stress testing mach 
####
#### doall loop
####
if test $DOALL -ge 1
then
    for FNAME in ${OFILE} ${OFILE}_all
    do
        echo "FNAME=${FNAME}"
        # clean up durn -9 that freaks MACH out as at july 2010
        sed 's/-9/0/g' $FNAME.ped > $FNAME.ped.clean
        rm -rf $FNAME.ped
        mv $FNAME.ped.clean $FNAME.ped 
        plink --file $FNAME --out $FNAME --nonfounders --freq --hardy2
        rm -rf $FNAME.dat
        echo "A Affection" > $FNAME.dat.tmp
        cut -f 2 $FNAME.map >> $FNAME.dat.tmp
        sed 's/rs/M rs/g ; s/RS/M RS/g' $FNAME.dat.tmp >> $FNAME.dat
        rm -rf $FNAME.dat.tmp
        MF=$FNAME
        RND=$RANDOM
        echo "seed=${RND}"
        if test $DEBUG -le 0
        then
            mach1 -d ${MF}.dat -p ${MF}.ped -s ${MF}.snps -h ${MF}.haps --greedy -r $ROU -o $MF --seed $RND >& ${MF}.stage1.log
        else
            mach1 -d ${MF}.dat -p ${MF}.ped -s ${MF}.snps -h ${MF}.haps --greedy -r $ROU -o $MF --seed $RND
        fi
        RND=$RANDOM
        echo "seed=${RND}"
        if test $DEBUG -le 0
        then
            mach1 -d ${MF}.dat -p ${MF}.ped -s ${MF}.snps -h ${MF}.haps --errorMap ${MF}.erate --crossoverMap ${MF}.rec --greedy --autoFlip --seed $RND --mle --mldetails -o $MF >& ${MF}.impute.log 
        else
            mach1 -d ${MF}.dat -p ${MF}.ped -s ${MF}.snps -h ${MF}.haps --errorMap ${MF}.erate --crossoverMap ${MF}.rec --greedy --autoFlip --seed $RND --mle --mldetails -o $MF 
        fi
    done
fi 
####
#### doall
####
# another pass
for FNAME in ${OFILE} ${OFILE}_all
do
    for NAME in FOUNDERS OFFSPRING
    do
        for FLAG in --filter-founders --filter-nonfounders
        do
            plink --file $FNAME $FLAG --recode --out ${FNAME}${NAME}
            sed 's/-9/0/g' ${FNAME}${NAME}.ped > ${FNAME}${NAME}.tmp
            rm -rf ${FNAME}${NAME}.ped
            mv ${FNAME}${NAME}.tmp ${FNAME}${NAME}.ped
            plink --file ${FNAME}${NAME} --out ${FNAME}${NAME} --nonfounders --freq --hardy2
            rm -rf ${FNAME}${NAME}.dat
            echo "A Affection" > ${FNAME}${NAME}.dat.tmp
            cut -f 2 ${FNAME}${NAME}.map >> ${FNAME}${NAME}.dat.tmp
            sed 's/rs/M rs/g ; s/RS/M RS/g' ${FNAME}${NAME}.dat.tmp >> ${FNAME}${NAME}.dat
            rm -rf ${FNAME}${NAME}.dat.tmp
            MF=${FNAME}${NAME}
            RND=$RAND
            echo "seed=${RND}"
            if test $DEBUG -le 0
            then
                mach1 -d ${MF}.dat -p ${MF}.ped -s ${FNAME}.snps -h ${FNAME}.haps --greedy -r $ROU -o $MF --seed $RND >& ${MF}.stage1.log
            else
                mach1 -d ${MF}.dat -p ${MF}.ped -s ${FNAME}.snps -h ${FNAME}.haps --greedy -r $ROU -o $MF --seed $RND 
            fi
            RND=$RANDOM
            echo "seed=${RND}"
            if test $DEBUG -le 0
            then
                mach1 -d ${MF}.dat -p ${MF}.ped -s ${FNAME}.snps -h ${FNAME}.haps --errorMap ${MF}.erate --crossoverMap ${MF}.rec --greedy --autoFlip --seed $RND --mle --mldetails -o $MF >& ${MF}.impute.log 
            else
                mach1 -d ${MF}.dat -p ${MF}.ped -s ${FNAME}.snps -h ${FNAME}.haps --errorMap ${MF}.erate --crossoverMap ${MF}.rec --greedy --autoFlip --seed $RND --mle --mldetails -o $MF
            fi
        done
    done
done



was

# experimental subset for mach family vs independent testing
# this selects 358 markers from chr22 care_sharp and then makes a random subset of 40% of the markers
#python makemachsamp.py /mnt/memefs/GWAS/SHARP/CARE_SHARP/CARE_SHARP testCARE 20000000 2000000 22 0.3
#for FNAME in testCARE testCARE_all
python makemachsamp.py testCARE_all test2 20000000 2000000 22 0.3
for FNAME in test2 test2_all
do
    plink --file $FNAME --out $FNAME --nonfounders --freq --hardy2
    rm -rf $FNAME.dat
    echo 'A Affection' > $FNAME.dat.tmp
    cut -f 2 $FNAME.map >> $FNAME.dat.tmp
    sed 's/rs/M rs/g ; s/RS/M RS/g' $FNAME.dat.tmp >> $FNAME.dat
    rm -rf $FNAME.dat.tmp
done
for FNAME in testCARE testCARE_all
do
    for NAME in FOUNDERS OFFSPRING
    do
        for FLAG in --filter-founders --filter-nonfounders
        do
            plink --file $FNAME --out ${FNAME}${NAME} $FLAG --recode
            plink --file ${FNAME}${NAME} --out ${FNAME}${NAME} --nonfounders --freq --hardy2
            rm -rf ${FNAME}${NAME}.dat
            echo 'A Affection' > ${FNAME}${NAME}.dat.tmp
            cut -f 2 ${FNAME}${NAME}.map >> ${FNAME}${NAME}.dat.tmp
            sed 's/rs/M rs/g ; s/RS/M RS/g' ${FNAME}${NAME}.dat.tmp >> ${FNAME}${NAME}.dat
            rm -rf ${FNAME}${NAME}.dat.tmp
        done
    done
done

"""
import sys,os,random,copy

debug = 1
g1000 = '/share/shared/data/1000g/2010_03' # 2010_06 seems to have errors?

def slicehaps(newPrefix=None,chrom='22',offset=20000000,length=2000000):
    """ take a slice from 1000g data corresponding to our slice
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
    snpf = '%s/snps/chr%s.snps' % (g1000,chrom)
    assert os.path.isfile(snpf),'# error - slicehaps cannot open snp file %s' % snpf
    hapf = '%s/hap/chr%s.hap' % (g1000,chrom)
    assert os.path.isfile(hapf),'# error - slicehaps cannot open hap file %s' % hapf
    mapf = '%s/map/chr%s.map' % (g1000,chrom)
    assert os.path.isfile(mapf),'# error - slicehaps cannot open map file %s' % mapf
    mapl = open(mapf,'r').readlines() # chr name offs
    mapls = [x.split('\t') for x in mapl]
    mapi = [i for i,x in enumerate(mapls) if (int(x[2]) >= offset and int(x[2]) <= (offset+length))] # use indices
    assert len(mapi) > 0, '# error slicehaps %s has zero length from %d to %d in %s' % (mapf,offset,offset+length,snpf)
    snps = open(snpf,'r').readlines()
    outsnps = [snps[i] for i in mapi] # subset	
    outsnpf = '%s.snps' % newPrefix
    outsnp = open(outsnpf,'w')
    outsnp.write(''.join(outsnps))
    outsnp.close()
    outhapf = '%s.haps' % newPrefix
    outhap = open(outhapf,'w')
    h = open(hapf,'r')
    for row in h:
        subh = [row[i] for i in mapi]
        subh.append('\n')
        outhap.write(''.join(subh)) 
    outhap.close()
    h.close()
    # now make links for _all files
    target = '%s_all' % newPrefix
    for e in ('.snps','.hap'):
        s = '%s%s' % (newPrefix,e)
        d = '%s%s' % (target,e)
        if not (os.path.isfile(d) or os.path.islink(d)):
            os.system('ln -s %s %s' % (s,d))

def slicePed(prefix=None,newPrefix=None,offset=0,length=100,chrom='22',randfrac=0.4,onlyOne=True):	
    """ take a slice from the map and pedfile
    produce an all subjects output ped file plus a complete trios output ped file
    at present multiple sibs will be output if onlyOne=False,
    because mach seems to be able to deal with these
    """
    assert prefix <> None,'slicePed error - no prefix for input'
    mfn = '%s.map' % prefix
    assert os.path.isfile(mfn)
    pfn = '%s.ped' % prefix
    assert os.path.isfile(pfn)
    mapl = open(mfn,'r').readlines()
    mapl = [x for x in mapl if len(x.split()) >=4]
    lm = len(mapl) # #snps
    mapls = [x.split() for x in mapl] # list of lists
    mslice = [i for i,x in enumerate(mapls) if (x[0]==chrom and int(x[3]) >= offset and int(x[3]) <= (offset+length))] 
    # keep this slice of markers
    assert len(mslice) > 1
    newmap = [mapl[x] for x in mslice] # haha
    cols = [6+2*x for x in mslice]
    cols += [7+2*x for x in mslice]
    cols.sort()
    if debug:
        print 'mslice[:20]=',mslice[:20]
    pfile = open(pfn,'r')
    newpfile = open('%s_all.ped' % newPrefix,'w')
    peds = {}
    offspring = {}
    trioped = {}
    noffs = {}
    affofs = {}
    fidDone = {}
    for row in pfile:
        srow = row.split()
        assert len(srow) == 6+2*lm,'slicePed error - short row len %d expected %d' % (len(srow),6+2*lm)
        ped = srow[:6]
        fid = ped[0]
        iid = ped[1]
        k = '%s_%s' % (fid,iid)
        peds.setdefault(k,copy.copy(ped))
        if ped[2] <> '0' and ped[3] <> '0':
            offspring.setdefault(k,copy.copy(ped))
            noffs.setdefault(fid,0)
            noffs[fid] += 1 # record each offspring
            if onlyOne:
                if ped[5] == '2': # affected
                    affofs.setdefault(fid,[])
                    affofs[fid].append(iid)
    pfile.seek(0)
    newmapf = open('%s_all.map' % newPrefix,'w')
    newmapf.write(''.join(newmap))
    newmapf.close()
    dotrio = False
    if len(offspring) > 0:
        dotrio = True
        triofile = open('%s_trios_all.ped' % newPrefix,'w')
    if len(offspring) > 0:
        for ped in offspring.values(): # flag offspring for inclusion or not
            fid,iid,pid,mid = ped[:4]
            dadk = '%s_%s' % (fid,pid)
            mumk = '%s_%s' % (fid,mid)
            kidk = '%s_%s' % (fid,iid)
            dad = peds.get(dadk,None)
            mum = peds.get(mumk,None)
            kid = peds.get(kidk,None)        
            if (dad != None) and (mum != None) and (kid != None):
                saveme = 1 # save all offspring if non onlyOne
                if onlyOne:
                    if fidDone.get(fid,None):
                        saveme = 0 # only one per family
                    else:
                        a = affofs.get(fid,None) # see if there are any
                        if a <> None:
                            if iid not in a: # there are affected and this is not one don't write
                                saveme = 0
                if saveme:
                    fidDone.setdefault(fid,1) # flag no more
                    trioped.setdefault(dadk,1) # include this complete trio in _trio file
                    trioped.setdefault(mumk,1)
                    trioped.setdefault(kidk,1)
    for row in pfile:
        srow = row.split()
        ped = srow[:6]
        k = '%s_%s' % (ped[0],ped[1])
        g = [srow[x] for x in cols]
        ped += g
        newpfile.write(' '.join(ped))
        newpfile.write('\n')
        if dotrio: 
            if trioped.get(k,None): # put out to triofile
                triofile.write(' '.join(ped))
                triofile.write('\n')
    pfile.seek(0)
    newpfile.close()
    if dotrio:
        triofile.close()
        newmapf = open('%s_trios_all.map' % newPrefix,'w')
        newmapf.write(''.join(newmap))
        newmapf.close()
    
    # now the random subset for testing
    n = int(len(mslice)*randfrac)
    assert n > 0,'slicePed error %f of %d < 1' % (randfrac,len(mslice))
    randm = random.sample(mslice,n) # subset of marker indices
    if debug:
        print 'presortrandm[:20]=',randm[:20]
    randm.sort() # sample should retain order but..doesn't seem to?
    if debug:
        print 'postsortrandm[:20]=',randm[:20]
    rcols = [6+2*x for x in randm] 
    rcols += [7+2*x for x in randm]
    rcols.sort() # yuck but gives the needed columns for each allele in order
    randmap = [mapl[x] for x in randm]
    newmapf = open('%s.map' % newPrefix,'w')
    newmapf.write(''.join(randmap))
    newmapf.close()
    if len(offspring) > 0:
        dotrio = True
        triofile = open('%s_trios.ped' % newPrefix,'w')
    newpfile = open('%s.ped' % newPrefix,'w')
    for row in pfile:
        srow = row.split()
        ped = srow[:6]
        k = '%s_%s' % (ped[0],ped[1])
        g = [srow[x] for x in rcols] # corresponding cols
        ped += g
        newpfile.write(' '.join(ped))
        newpfile.write('\n')
        if dotrio:
            if trioped.get(k,None): # put out to triofile
                triofile.write(' '.join(ped))
                triofile.write('\n')
    newpfile.close()
    pfile.close()
    if dotrio:
        triofile.close()
        newmapf = open('%s_trios.map' % newPrefix,'w')
        newmapf.write(''.join(randmap))
        newmapf.close()
    if debug:
        print 'slicing %s from %d to %d' % (pfn,cols[0], cols[-1])
        print 'random cols = %s' % rcols
    
    
"""for FNAME in testCARE testCARE_all
do
    for NAME in FOUNDERS OFFSPRING
    do
        for FLAG in --filter-founders --filter-nonfounders
        do
            plink --file $FNAME --out ${FNAME}${NAME} $FLAG --recode
            plink --file ${FNAME}${NAME} --out ${FNAME}${NAME} --nonfounders --freq --hardy2
            rm -rf ${FNAME}${NAME}.dat
            echo 'A Affection' > ${FNAME}${NAME}.dat.tmp
            cut -f 2 ${FNAME}${NAME}.map >> ${FNAME}${NAME}.dat.tmp
            sed 's/rs/M rs/g ; s/RS/M RS/g' ${FNAME}${NAME}.dat.tmp >> ${FNAME}${NAME}.dat
            rm -rf ${FNAME}${NAME}.dat.tmp
        done
    done
done
"""


class machRun():
    """
    mach run container - inputs, intermediates, outputs, settings 
    """
    
    
    def __init__(self, basename=None, refname=None, stage1name=None, outname=None,rounds=None,inseed=0,machexe='mach1',**kwd):
        """
        basename is ped file basename
        refname is basename for the set of reference haplotypes - hap/snps pair
        stage1name is the name for the crossover and erate files from mach stage 1 - made if not exists
        outname is used as stage2 output name
        """
        self.machexe = machexe
        self.basename = basename     
        self.refname = refname
        self.rounds = rounds
        self.inseed = seed
        self.stage1name = stage1name
        self.stage2name = outname
        self.crossratename = stage1name
        self.eratename = stage1name
        assert rounds > 0,'## error machRun class instance initiated with rounds not > 0'
        assert self.basename <> None,'## error machRun class instance initiated with missing basename'
        assert os.path.exists('%s.ped' % self.basename),'## error no ped file found for machRun class instance basename %s' % basename
        assert self.refname <> None,'## error machRun class instance initiated with missing reference haplotype name'
        assert os.path.exists('%s.snps' % self.refname),'## error no snps file found for machRun class instance refname %s' % refname
        assert os.path.exists('%s.hap' % self.refname),'## error no hap file found for machRun class instance refname %s' % refname
        datname = '%s.dat' % self.basename
        if not os.path.exists(datname):
            self.makedat()
        if not (os.path.exists('%s.erate' % stage1name) and os.path.exists('%s.rec' % stage1name)):
            doMach1(self)
        doMach2(self)
        
    def nextSeed(self):
        if self.inseed <> 0:
            self.seed = self.inseed
        else:
            self.seed = random.randint(1,sys.maxint)
            
    def makeDat(self):
        """ make a dumb old .dat file for mach to match the map - note 
        additional affection pheno field so mach can read plink style 6 column pedigrees
        """
        mapname = '%s.map' % self.basename
        assert os.path.isfile(mapname) or os.path.islink(mapname)
        m = open(mapname,'r').readlines()
        mlist = [x.strip().split() for x in m if len(x.strip().split()) >=3]
        rslist = ['M %s' % x[1] for x in mlist]
        datname = '%s.dat' % self.basename
        os.system('rm -rf %s' % datname)
        outf = file(datname,'w')
        outf.write('A Affection\n')
        outf.write('\n'.join(rslist))
        outf.write('\n')
        outf.close()        
        
    def doMach1(self):
        """
        stage1 - given a few hundred independent subjects for a whole or part chromosome, corresponding regional snp and hap data 
        make crossovermap and error rates for stage 2
        
        expects a ped file at basepath

        mach1 -d ${MF}.dat -p ${MF}.ped -s ${FNAME}.snps -h ${FNAME}.haps -r $ROU -o $MF --seed $RND >& ${MF}.stage1.log
        """
        self.seed = self.nextseed()
        if not os.path.exists('%s.dat' % self.basename):
            self.makedat()
        cl = 'mach1 -d %s.dat -p %s.ped -s %s.snps -h %s.haps -r %d -o %s -seed %d >& %s_mach1.log' % (self.basename,self.basename,
          self.refname,self.refname,self.rounds,self.stage1name,self.seed,self.stage1name) 
        print '## doMach1 cl = ',cl
        
        
    def doMach2(self):
        """
        need basename=None,refname=None, eratename=None, crossratename=None, outname=None,haps=None,snps=None,rounds=None,seed=0
        stage2
        mach1 -d ${MF}.dat -p ${MF}.ped -s ${FNAME}.snps -h ${FNAME}.haps --errorMap ${MF}.erate 
        --crossoverMap ${MF}.rec --greedy --autoFlip --seed $RND --mle --mldetails -o $MF >& ${MF}.impute.log 
        if no erate/crossrate files exist, make them - this could take a long time
        """
        if not (os.path.exists(self.eratename) and os.path.exists(self.crossratename)): # make
            self.seed = self.nextseed()
            doMach1(self)
        self.seed = self.nextseed()
        cl = """mach1 -d %s.dat -p %s.ped -s %s.snps -h %s.haps --errorMap %s.erate --crossoverMap %s.rec --greedy --autoFlip --seed %d --mle --mldetails
         -o %s >& %s_mach2.log""" % (self.basename,self.basename,self.refname,self.refname,self.stage1name,self.stage1name,self.seed,self.outname, self.outname)
        print '## doMach2 cl = ',cl

class plinkRun():
    """
    plink run container - inputs, intermediates, outputs, settings 
    """   
    def __init__(self, **kwd):
        self.plinkPath = kwd.get('plinkPath','plink')
        self.baseName=kwd.get('baseName',None)
        self.baseFType=kwd.get('baseFType',None)
        self.extraParams = kwd.get('extraParams',None)
        self.outPath=kwd.get('outPath',None)
        self.logPath=kwd.get('outLogPath','plink.log')
        self.plinkTasks=kwd.get('plinkTasks',[])
        self.vclBase=kwd.get('vclBase','') # whatever it takes at the start of a cl to invoke plink - often 'plink'
        self.cd=kwd.get('cd','')
        self.clBase=kwd.get('clBase','')
        self.setBase = None
        if self.baseName:
            if self.baseFType:
                if self.baseFType == 'pbed':
                    self.setBase = '--bfile %s' % self.baseName
                else:
                    self.setBase = '--file %s' % self.baseName
        if self.outPath:
            self.setOut = '--out %s' % self.outPath        
            
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
       
    def addFilter(self,ftype='--filter-founders',extraOptions=''):
        assert ftype in ['--filter-founders','--filter-offspring'],'#error plinkSettings addfilter - ftype = %s' % ftype
        if len(extraOptions) > 0:
            extraOptions = '%s %s' % (extraOptions,ftype)
        else:
            extraOptions = ftype
        self.doThese(extraOptions=extraOptions)
        
    def run(self):
        """
        plink blathers when doing pruning - ignore
        
        """
        fplog,plog = tempfile.mkstemp()
        alog = []
        alog.append('## Rgenetics: http://rgenetics.org Galaxy Tools rgQC.py Plink runner\n')
        for task in self.plinkTasks: # each is a list
            vcl = self.vclbase + task
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
        outname = names[i]
        p = plinkRun(baseName=basename,baseFType='lped',outName=outname)
        p.addFilter(ftype=flag,extraOptions='')
        p.doThese(extraOptions='--nonfounders --freq --hardy2')
        alog = p.run()

def machtest():
    npar = 9
    if len(sys.argv) < npar:
        print 'Usage: makemachsamp.py infile_prefix outfile_prefix offset slicelength chrom randfrac onlyOne inseed'
        sys.exit(1)
    prefix,outprefix,offset,length,chrom,randfrac,onlyOne,inseed = sys.argv[1:npar]
    onlyOne = onlyOne and onlyOne.lower() in ['true','1']
    if inseed == '':
        inseed = 0
    else:
        try:
            inseed = int(inseed)
        except:
            print 'inseed %s not an integer in cl = %s' % (inseed,sys.argv[1:npar])
    slicePed(prefix=prefix,newPrefix=outprefix,offset=int(offset),length=int(length),chrom=chrom,randfrac=float(randfrac),onlyOne=True)
    slicehaps(newPrefix=outprefix,chrom=chrom,offset=int(offset),length=int(length))
    makeFFO(basename=outprefix) # create offspring/founders only ped files for mach phase 1
    fped = '%sFOUNDERS' % outprefix
    m = machRun(basename=fped,refname=outprefix,stage1name=outprefix,outname=outprefix,rounds=int(rounds),seed=inseed,machexe='mach1')
    

if __name__ == "__main__":
    machtest()



