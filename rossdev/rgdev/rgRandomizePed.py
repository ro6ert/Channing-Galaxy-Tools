"""
beautiful but useless
In caco, you only have allele frequency differences, and
allele sharing gives more information about relatedness.
Allele sharing in segments?
#
#
but only at margins
added filter for autosomes - getting 10^-11 p values for icondb unfiltered!
added ALL flag to convert families into all unrelateds
# for R - for stoopid R to read plink assoc files, MUST use whitespace option 
# eesh.
a = read.delim(file='test.assoc',header=T,strip.white=T,sep="")  
pvec = a$P
npvals = length(a$P)
qqplot(-log10(pvec),-log10(runif(npvals)))      
mx = c(0,log10(npvals))        
points(mx,mx,type='l')  
# a = read.delim(file='test.assoc',header=T,strip.white=T,sep="")
# randomize genotypes - eg icondb
# permute subjects for each locus so allele frequencies are
# unchanged - ld will be stuffed but too bad..
# remember to shuffle affection and change famid
# data comes with snp as cols
# so need a row-wise generator - random reads to file?
"""

import random,sys,copy,shutil,time
from scipy import weave
from numpy import array


debug = 1 # turn off for galaxy - prints do not show

# thinking about weave, but numpy will probably do
ccode="""
/*
        for snp in self.useme: # sample colwise to create reverse affection subjects
            if (snp+1) % 50000 == 0:
                dur = time.time() - started
                logf.write('randomizing at snp %d of %d; %d alleles/sec\n' % (snp,
                                    self.nsnp,int(snp*self.nped/dur)))
            allele1 = 2*snp # offset
            randvec = []
            for sampi in range(self.nped): # for each sample
                samp = self.genos[sampi] # get genos
                randsamp = random.choice(gchoice) # choose another subject at random
                x = self.genos[randsamp] # choose a subject at random
                randvec.append([copy.copy(x[allele1]),copy.copy(x[allele1 + 1])]) # randomized value for this persons genotype
            for sampi in range(self.nped): # now have a vector of random genos sampled from identical allele freq dist
                samp = self.genos[sampi]
                samp[allele1] = randvec[sampi][0]
                samp[allele1 + 1] = randvec[sampi][1]

need useme,nped,nsnp,genos passed in
*/
#include "stdlib.h"

char randvec[nped*2];
unsigned int now = time(NULL);
srand( now ); /* seed rng */
for (snpi=0; snpi<nsnp; snpi++) {
    int snp = useme[snpi];
    int allele1 = 2*snp;
    for (sampi=0; sampi<nped; sampi++) {
        int randsamp = rand() % nped; /* random with replacement subject */
        randvec[2*sampi] = genos[randsamp][allele1];
        randvec[2*sampi + 1] = genos[randsamp][allele1+1];
        }
    for (sampi=0; sampi<nped; sampi++) {
        genos[sampi][allele1] = randvec[2*sampi];
        genos[sampi][allele1+1] = randvec[2*sampi+1];        
        }
    }
return_val = genos;
"""

class randped:
    """ tricky class to write a caco data set
    adding a new record with opposite affection status for every subject
    and randomized genotypes by locus so allele freqs are identical
    to seek bias from relatedness in caco studies - eg illumina icondb
    """
    def __init__(self,pedf=None,mapf=None,logf=None,useAll=0):
        self.logf = logf
        genos = []
        pedigrees = []
        pedigreesrev = []
        autosomes = map(str,range(1,23))
        autodict = dict(zip(autosomes,autosomes)) # quick check
        self.maplist = mapf.readlines() # for writing later
        maplist = [x.split() for x in self.maplist if (len(x.split()) > 0)]
        chrlist = [x[0] for x in maplist]
        useme = [i for i,x in enumerate(chrlist) if autodict.get(x,None)] # autosomal markers
        self.nsnp = len(useme)
        self.nalleles = 2*self.nsnp
        usealleles = [] # need actual allele places for speed
        for i in useme:
            usealleles.append(2*i)
            usealleles.append(2*i + 1) # we need alleles
        self.useme = useme
        started = time.time()
        fakeFID = {}
        if debug:
            self.logf.write('using %d autosomal markers of %d read\n' % (len(self.useme),len(maplist)))
        for n,row in enumerate(pedf):
            if debug and (n+1) % 1000 == 0:
                dur = time.time() - started
                self.logf.write('processing row %d - %f secs per row\n' % (n,dur/n)) 
            lrow = row.strip().split()
            if len(lrow) >= 8:
                founder = (lrow[2] == '0' and lrow[3] == '0') # founders only usually
                if founder or useAll:
                    p = lrow[:6]
                    if useAll and not founder: # fake family id to make 'unrelated'
                        p[0] = '%s_%d' % (p[0],n)
                        p[2] = '0'
                        p[3] = '0'
                    p2 = copy.copy(p)
                    g = lrow[6:]
                    guseme = ['%s %s' % (g[2*i],g[2*i +1]) for i in useme] # only autosomal marker alleles
                    # we need to shuffle genotype pairs and split them for printing
                    #grow = array.array('c',guseme)
                    grow = array(guseme) # numpy more better
                    genos.append(grow) # rest as an array to save storage
                    pedigrees.append(copy.copy(p)) # pedigree
                    if p2[5] == '1':
                        p2[5] = '2' # reverse affection for shuffled genome
                    else:
                        p2[5] = '1'
                    pedigreesrev.append(p2) # reversed affection
        self.nped = len(pedigrees)
        self.logf.write('found %d snp and %d subjects\n' % (self.nsnp,self.nped))
        self.genos = genos
        self.pedigrees = pedigrees
        self.pedigreesrev = pedigreesrev
                

    def slownextRow(self):
        """
        """
        for n,g in enumerate(self.genos): # for each subject's geno
            row = copy.copy(self.pedigrees[n])
            row += ''
            row += list(g)
            yield '\t'.join(row)        
        gchoice = xrange(self.nped)
        started = time.time()
        for snp in self.useme: # sample colwise to create reverse affection subjects
            if (snp+1) % 50000 == 0:
                dur = time.time() - started
                self.logf.write('randomizing at snp %d of %d; %d alleles/sec\n' % (snp,
                                    self.nsnp,int(snp*self.nped/dur)))
            allele1 = 2*snp # offset
            randvec = []
            for sampi in range(self.nped): # for each sample
                samp = self.genos[sampi] # get genos
                randsamp = random.choice(gchoice) # choose another subject at random
                x = self.genos[randsamp] # choose a subject at random
                randvec.append([copy.copy(x[allele1]),copy.copy(x[allele1 + 1])]) # randomized value for this persons genotype
            for sampi in range(self.nped): # now have a vector of random genos sampled from identical allele freq dist
                samp = self.genos[sampi]
                samp[allele1] = randvec[sampi][0]
                samp[allele1 + 1] = randvec[sampi][1]
        for n,g in enumerate(self.genos): # for each random/reversed subject's geno
            row = copy.copy(self.pedigreesrev[n])
            row[0] = '%d%s' % (n+1,row[0]) # make fake fid
            row[1] = '%d%s' % (n+1,row[1]) # make fake iid
            row += ''
            row += g
            yield '\t'.join(row)        

    def nextRow(self):
        """
        numpy is much better
>>> foo = [array(['%d' % random.randint(1,9) for x in range(10)]) for y in range(5)]
>>> bar = array(foo)
>>> bar
array([['8', '8', '7', '2', '5', '4', '7', '7', '6', '1'],
       ['2', '8', '4', '9', '2', '6', '4', '3', '2', '9'],
       ['6', '2', '1', '7', '3', '6', '3', '5', '4', '8'],
       ['3', '3', '6', '9', '6', '5', '4', '6', '4', '7'],
       ['9', '2', '9', '6', '9', '3', '4', '5', '9', '1']], 
      dtype='|S1')
>>> random.shuffle(bar[:,0])
>>> bar
array([['6', '8', '7', '2', '5', '4', '7', '7', '6', '1'],
       ['8', '8', '4', '9', '2', '6', '4', '3', '2', '9'],
       ['2', '2', '1', '7', '3', '6', '3', '5', '4', '8'],
       ['9', '3', '6', '9', '6', '5', '4', '6', '4', '7'],
       ['3', '2', '9', '6', '9', '3', '4', '5', '9', '1']], 
      dtype='|S1')
>>> random.shuffle(bar[:,2])
>>> bar
array([['6', '8', '1', '2', '5', '4', '7', '7', '6', '1'],
       ['8', '8', '7', '9', '2', '6', '4', '3', '2', '9'],
       ['2', '2', '6', '7', '3', '6', '3', '5', '4', '8'],
       ['9', '3', '9', '9', '6', '5', '4', '6', '4', '7'],
       ['3', '2', '4', '6', '9', '3', '4', '5', '9', '1']], 
      dtype='|S1')
        """
        for n,g in enumerate(self.genos): # for each subject's geno
            row = copy.copy(self.pedigrees[n])
            row += ''
            row += list(g)
            yield ' '.join(row) # g is already space delim allele pairs        
        started = time.time()
        genos = array(self.genos) # convert to numeric
        for snp in range(self.nsnp): # sample colwise to create reverse affection subjects
            if (snp+1) % 50000 == 0:
                dur = time.time() - started
                self.logf.write('randomizing at snp %d of %d; %d alleles/sec\n' % (snp,
                                    self.nsnp,int(snp*self.nped/dur)))
            gsamp = genos[:,snp] # vector of sample genos at this locus
            rv = []
            for ped in range(self.nped): # now need to randomly sample and replace
                rsamp = random.randint(0,self.nped-1) # any subject's genotype pair
                rv.append(gsamp[rsamp]) # grab a subject with replacement
            genos[:,snp] = array(rv) # replace slice with a randomized version
        for n,g in enumerate(genos): # for each random/reversed subject's geno
            row = copy.copy(self.pedigreesrev[n])
            row[0] = '%d%s' % (n+1,row[0]) # make fake fid
            row[1] = '%d%s' % (n+1,row[1]) # make fake iid
            row += ''
            row += g
            yield ' '.join(row)        
        

    def write(self,f=None):
        """
        """
        for row in self.nextRow():
            f.write(row)
            f.write('\n')
            
    def writemap(self,mapf=None):
        """
        """
        f = file(mapf,'w')
        for i in self.useme: # indexes
            f.write('%s\n' % self.maplist[i])
        f.close()            
            
def test():
    """
    borked as write won't work but used to develop
    """
    sp = """1 1 0 0 1 2 1 2 3 4\n1 2 0 0 2 1 1 2 3 4\n1 3 1 2 1 2 2 2 4 4\n"""
    p = sp.split('\n')
    q = []
    foo = randped(p,None)
    foo.write(q)
    for row in q:
        print row
    
def rPed():
    """
    randomize - call with pedpath outpath
    """
    v = sys.argv
    if len(v) >= 1:
        pedname = v[1]
        pedpath = '%s.ped' % pedname
        mappath = '%s.map' % pedname
        logpath = '%s_randomised.log' % pedname
        outpath = '%s_randomised.ped' % pedname
        outmap = '%s_randomised.map' % pedname
        outfile = file(outpath,'w')
        useAll = 0
        if len(v) >= 2:
            if v[1].upper() == 'ALL':
                useAll = 1
                print 'Pretending families are unrelated'
        p = file(pedpath,'r')
        m = file(mappath,'r')
        logf = file(logpath,'w')
        rp = randped(p,m,logf)
        rp.write(outfile)
        outfile.close()
        rp.writemap(outmap)
        logf.close()
    else:
        print 'need pedname - will randomize pedname.ped to pedname_randomized.ped'
        
if __name__ == "__main__":
    rPed()
