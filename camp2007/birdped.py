# birdped.py
# fix affy 6.0 sample birdseed calls with fake pedigree info
# the calls are -9,0,1,2 for number of ?rare alleles

# eesh - duplicate rs numbers are output in the genotype file!

import csv,sys,os,array

samptable = 'NHLBI_Sample_Table.txt'
genofile = 'birdseed.calls.txt'
annotable = 'GenomeWideSNP_6.na24.annot.csv'
alleletran = {'A':'1','C':'2','G':'3','T':'4','N':'0'}

def getSamp(stf=samptable):
    """
    Plate       Well    Plate:Well      Sample Name     Best Array      QC cr   birdseed cr     Trio Info       Trio Accuracy   Pico Green Conc. 
    SMP4_0002074        A01     SMP4_0002074:A01        CA001557609DD   a520532-00-651208-082008-4034931-60498  93.98000336     98.061  1037_Father
97.514  48.2608527
    SMP4_0002074        A02     SMP4_0002074:A02        CA001558680DD   a520532-00-651208-082008-4034931-60494  96.62000275     98.05   1037_Mother
97.514  48.74598601

    """
    f = open(stf,'r')
    st = f.readlines()
    sl = [x.strip().split('\t') for x in st]
    head = [x.replace(' ','_') for x in sl[0]]
    sl = [x for x in sl[1:] if len(x) == len(head)] # get rid of duff stuff at end
    caid = head.index('Best_Array')
    ctrio = head.index('Trio_Info')
    aids = [x[caid] for x in sl] # list of array ids
    trios = [x[ctrio] for x in sl] # list of stuff about sample
    peds = []
    roles = ['Father','Mother','Proband']
    for t in trios: # create fake peds
        status = t.split('_')
        if len(status) >= 2:
           x = id,role = t.split('_')[:2]
           fid = id
           iid = roles.index(role) + 1 # fake
           if iid == 3:
              s = '%s %d 1 2 1 2 ' % (fid,iid)
           else:
              s = '%s %d 0 0 %d 1' % (fid,iid,iid) # doubles as gender for parents only
        else:
            s = 'hapmap'
        peds.append(s)
    pdict = dict(zip(aids,peds)) # lookup
    print 'pdict=',pdict
    return pdict

def getAnno(annf=annotable):
    """
 "Probe Set ID","Affy SNP ID","dbSNP RS ID","Chromosome","Physical Position","Strand","ChrX pseudo-autosomal region 1",
 "Cytoband","Flank","Allele A","Allele B","Associated Gene","Genetic Map","Microsatellite","Fragment Enzyme Length Start Stop",
 "Allele Frequencies","Heterozygous Allele Frequencies","Number of individuals/Number of chromosomes",
 "In Hapmap","Strand Versus dbSNP","Copy Number Variation","Probe Count","ChrX pseudo-autosomal region 2","In Final List",
 "Minor Allele","Minor Allele Frequency"   

    we make a dynamic dict so 9,0,1,2 can be translated into a genotype
    """
    dictanno = {}
    dictallele = {}
    rsdict = {} # so we can renumber duplicate rs numbers - otherwise rgQC gets unhappy
    f = csv.reader(open(annf,'r'))
    dups = []
    for i,row in enumerate(f):
        if row[0][0] <> '#' and row[0] <> 'Probe Set ID': # not header or comment
            affyid = row[0]
            rs = row[2]
            if rs == '---':
                rs = affyid
            rsrepeat = rsdict.get(rs,0)
            if rsrepeat:
                rsdict[rs] += 1
                rs = '%s_%d' % (rs,rsrepeat)
                dups.append(rs) # for posterity
            else:
                rsdict[rs] = 1 # first time seen - prepare for the worst
            chrom = row[3]
            offset = row[4]
            strand = row[5]
            a = row[9]
            fbata = alleletran.get(a,'0')
            b = row[10]
            fbatb = alleletran.get(b,'0')
            gmap = row[12].split('//')[0].strip()
            dictanno[affyid] = [rs,chrom,offset,strand,a,b,gmap]
            # note we translate -1 into 9 so we can use arrays to store genotypes
            dictallele[affyid] = {'9':'0 0','0':'%s %s' % (fbata,fbata),
                                  '1':'%s %s' % (fbata,fbatb),'2':'%s %s' % (fbatb,fbatb)}
            if (i+1) % 100000 == 0:
              print '#getAnno at row %d in %s' % (i,annf) 
              print 'a=%s,b=%s, dictallele = %s for affyid %s' % (a,b,dictallele[affyid],affyid)
    return dictanno,dictallele,dups 
    

def transposeIt(dat=[]):
    """ generic transpose iterator - yields transposed rows of a rectangular list of lists of genotypes
    """
    nrec = len(dat) # number of markers
    if nrec == 0: # nothing to do
       yield []
       raise StopIteration
    else:
       dlen = len(dat[0]) # number of subjects*2 since each row has 2 alleles for each subject
       if dlen % 2 <> 0:
           print 'Odd transpose data row length %d - cannot extract pairs of alleles' % dlen
           yield []
           raise StopIteration
       for i,row in enumerate(dat): # check rectangularity
           if len(row) <> dlen:
               print '## transposeIt got an irregular data structure - row #%d length %d instead of %d' % (i+1,len(row),dlen)
               yield []
               raise StopIteration      
       print '## transposeIt got %d subjects, each with %d markers' % (dlen,nrec)
       for subject in range(dlen):
           res = array.array('c',['0',]*nrec) # empty array for each subject with a slot for every geno    
           for marker,d in enumerate(dat): # grab each marker in turn
               res[marker] = d[subject]
           yield res

def repDups(duprslist=[],gt=[],idlist=[],rsvec=[]):
    """some rs appear on both sty and nsp chips so duplicates - flagged as rsxxx_2
    """
    rsdict = dict(zip(rsvec,range(len(rsvec)))) # lookup rs to get row in untransposed gt 
    subjdisc = [0 for x in idlist]
    rsdisc = [0 for x in duprslist]
    diff = {}
    for n,rs in enumerate(duprslist):
      firstpos = rsdict.get(rs,None)
      realrs = '%s' % rs.split('_')[0]
      secondpos = rsdict.get(realrs,None)
      if not firstpos or not secondpos:
          print '### repDups problem at dup %d - realrs %s or rs %s not found in rsdict' % (n,realrs,rs)
      else:
          gt1 = gt[firstpos]
          gt2 = gt[secondpos]
          diff[rs] = [(i,x,gt2[i]) for i,x in enumerate(gt1) if x <> gt2[i]] 
          for (i,x,y) in diff[rs]:
              subjdisc[i] += 1
          rsdisc[n] += len(diff[rs])
    duprep = file('birdpedDups.txt','w')
    subjdiscrep = [(subjdisc[i],idlist[i]) for i,id in enumerate(idlist) if subjdisc[i] > 0]
    # ('1233',33)
    subjdiscrep.sort()
    subjdiscrep.reverse()
    subjdiscrep = ['%d discordances for subject %s' % x for x in subjdiscrep]
    s = '\n'.join(subjdiscrep)
    print s
    duprep.write(s)
    duprep.write('\n')
    rsdiscrep = [(rsdisc[i],duprslist[i]) for i,rs in enumerate(duprslist) if rsdisc[i] > 0]
    # ('rs123',27)
    rsdiscrep.sort()
    rsdiscrep.reverse()
    rsdiscrep = ['%d discordances for marker %s' % x for x in rsdiscrep]
    s = '\n'.join(rsdiscrep)
    print s
    duprep.write(s)
    duprep.write('\n')
    duprep.close()


def getGeno(gtf=genofile):
    """
probeset_id     cels/a520532-00-650360-082008-4034757-32470.CEL cels/a520532-00-650371-082008-4034757-32602.CEL 
cels/a520532-00-651208-082008-4034931-L
SNP_A-2131660   1       2       1       2       1       2       2       2       2       2       2       2       2       2       1       1       
2     2
SNP_A-1967418   2       2       2       2       1       2       2       1       2       2       2       1       2       2       1       1       
2     2
S
    """
    pdict = getSamp()
    annodict,alleledict,dupes = getAnno()
    if len(dupes) > 0:
        f = file('affydupes.xls', 'w')
        f.write('\n'.join(dupes))
        f.close()
        print '## wrote %d duplicate rs' % (len(dupes))
    f = open(gtf,'r')
    gh = None
    gl = []
    snpvec = []
    acolvec = [] # autosome column vector
    for i,row in enumerate(f):
        if row[0] <> '#':
            lrow = row.strip().split()
            if lrow[0] == 'probeset_id':
                hids = lrow[1:]
            else:
                snpid = lrow[0]
                snpvec.append(snpid)
                lrow = lrow[1:]
                lrow = [x.replace('-1','9') for x in lrow]
                newg = array.array('c',lrow)
                gl.append(newg)                
    pmapf = file('birdlped.map','w')
    amapf = file('birdlpedAutos.map','w')
    for n,affyid in enumerate(snpvec):
        rsa = annodict[affyid] # dictanno[affyid] = (rs,chrom,offset,strand,a,b,gmap)
        row = [rsa[1],rsa[0],rsa[6],rsa[2]]
        pmapf.write('\t'.join(row))
        pmapf.write('\n')
        try:
            c = int(rsa[1]) # autosomal chrom
            amapf.write('\t'.join(row))
            amapf.write('\n')
            acolvec.append(n) # write this geno to autosomal file
        except:
            pass
    pmapf.close()
    hids = [x.split('/')[1] for x in hids] # get rid of leading directory
    hids = [x.split('.')[0] for x in hids] # clean up .cel
    pedvec = [pdict.get(x,[]) for x in hids]    
    rsvec = [annodict.get(x,['unknownRS',])[0] for x in snpvec] # convert affysnpid to rs
    arsvec = [rsvec[x] for x in acolvec] # autosomal rs
    repDups(duprslist=dupes,gt=gl,idlist=pedvec,rsvec=rsvec)
    gtranvec = [alleledict.get(x,None) for x in snpvec] # dicts for converting to fbat alleles
    t = transposeIt(gl) # iterator to return transposed lists
    foutf = file('birdfbat.ped','w')
    afoutf = file('birdfbatAutos.ped','w')
    foutf.write(' '.join(rsvec))
    foutf.write('\n')
    afoutf.write(' '.join(arsvec)) # autosomes
    afoutf.write('\n')
    poutf = file('birdlped.ped','w')
    apoutf = file('birdlpedAutos.ped','w')
    for i,row in enumerate(t): # each row is the geno call (0,1,2 or 9 for missing)
        ped = pedvec[i]
        if ped <> 'hapmap':
            print '## at %d, have ped = %s' % (i,ped)
            # (0,1,2 or 9)
            geno = [gtranvec[n].get(x,'0 0') for n,x in enumerate(row)]
            if i < 5:
                print 'row = %s' % (' '.join(row[500:505]))
                print 'geno=',geno[1000:1100]
            # whew. Get the right alleles for each snp
            s = '%s %s\n' % (ped,' '.join(geno))
            foutf.write(s)
            poutf.write(s)
            geno = [gtranvec[n].get(row[n],'0 0') for n in acolvec] # autosomes only
            if i < 5:
                print 'autogeno=',geno[1000:1100]
            # whew. Get the right alleles for each snp
            s = '%s %s\n' % (ped,' '.join(geno))
            afoutf.write(s)
            apoutf.write(s)
        else:
            print 'hapmap sample %d ignored' % i
    foutf.close() 
    poutf.close()
    afoutf.close() 
    apoutf.close()

if __name__ == "__main__":
    getGeno()
