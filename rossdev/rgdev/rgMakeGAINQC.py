"""
convert framingham share to gainqc format
See http://www.sph.umich.edu/csg/abecasis/GainQC/index.html

share matrix format has 
[rerla@hg gainqc]$ head General_Release.geno.chr1  | more
dbSNP_ss dbSNP_rs Affy_SNPID chr pos ss2rs_orient rs2genome_orient allele1 allele2 11359 
"""
import sys,os,operator

missval = 'N'

AFFECTED = 'STATUS=AFFECTED'
UNAFFECTED = 'STATUS=UNAFFECTED'

def getPed(inped='phs000007.pht000183.v1.p1.c1.share_ped.GRU.txt'):
    """read the pedigree - framingham in this current instance
    """
    ped = {}
    pf = file(inped,'r')
    for n,row in enumerate(pf):
        if row[0] <> '#' and row[:5] <> 'pedno':
           lrow = row.split('\t')
           if len(lrow) >= 6:
              famId,shareId,paId,maId,gender = lrow[:5]
              if paId == '':
                  paId = '0'
              if maId == '':
                  maId = '0'
              if maId and paId:
                  affstatus = AFFECTED
              else:
                  affstatus = UNAFFECTED
              ped[shareId] = (famId,shareId,paId,maId,gender,affstatus)
           else:
              print 'dud row %s in pedfile' % lrow
    return ped

def getShareGenos():
    """ this one works for framingham share matrix format
    ross lazarus march 2008
    """
    for chromn,chrom in enumerate(chromlist):
         infiles = []
         heads = []
         sids = []
         snpAnno = {}
         genos = {}
         fnames = []
         for use in ['General_Release.geno.','NonProfit_Release.geno.']:
            fname = '%schr%s' % (use,chrom)
            fnames.append(fname)
            f = file(fname,'r')               
            head = f.readline().strip().split()
            heads.append(head)
            sid = head[9:]
            infiles.append(f)
            if chromn == 0: # stow firsttime, check after
               sampleIds += sid
            else:
               sids += sid # for checking
         if chromn <> 0: # are the samples in the same exact order?
            for j,id in enumerate(sids):
                if sampleIds[j] <> id:
                   print 'Matrix file %s has out of order sample ids!' % fname
                   sys.exit(1)                      
         for fnum,f in enumerate(infiles):
            f.seek(0)
            for i,row in enumerate(f): # rewind to avoid next/readline problems!
                if i > 0:
                    if i % 5000 == 0:
                          print 'At line %d in %s' % (i,fnames[fnum])
                    lrow = row.strip().split()
                    #dbSNP_ss dbSNP_rs Affy_SNPID chr pos ss2rs_orient rs2genome_orient allele1 allele2
                    rs = lrow[1]
                    chrom = lrow[3]
                    try:
                            pos = int(lrow[4])
                    except:
                            pos = 0
                    a1 = lrow[7]
                    a2 = lrow[8]
                    if not snpAnno.get(rs,None):
                         snpAnno[rs] = [rs,chrom,pos,'A','[%s/%s]' % (a1,a2)]     
                         #gainqc wants PREFERRED_ID   CHR   POSITION        QC_TYPE SEQUENCE
                         #             Marker1        1     3019921         A       CGGCT[A/T]ACGTA
                    g = genos.get(rs,[])
                    thisrow = lrow[9:]
                    for gtn,gt in enumerate(thisrow):
                         if gt == 'ND':
                            thisrow[gtn] = 'NN' # ugh. details.                                          
                    g += thisrow # extend - subjects are known same order
                    genos[rs] = g # append subjects to each row
    return genos, snpAnno



def getLpedGenos(basename='',repository='./'):
    """ for plink style linkage format map/ped files
    ross lazarus march 2008
    """
    
    def gainQCRowGenerator(destname=None,outclass=None):
        """
        for ped files
        eig is a row=marker, column=subject transposed format like gainqc wants
        except of course, different. Why can't we all just get along?
        This is brainfuck not python. Generates single marker
        genotype rows selected from among classes of subjects as specified in
        a list as outclass. The eig dict has each class's subjects genotypes
        all in the same marker order, so we can take the same row for any given marker

        """
        for n,row in enumerate(self.eig[outclass[0]]): # have to write fully row wise
            s = ''.join(row)
            if len(outclass) > 1: # must append cases to rows snpwise
                for aclass in outclass[1:]: # rest
                    appendme = ''.join(self.eig[aclass][n]) # take the nth row
                    s = '%s%s' % (s,appendme)
            yield s # we're a generator so iterable. Neat.


    pedroot = os.path.join(sourcedir,basename)
    pedf = '%s.ped' % pedroot
    f = file(pedf,'r')
    if linkage: # read map file
        map = readMap(pedroot)
        rslist = [x[1] for x in map] # get rs numbers
        outmap = file(mappath,'w')
        maps = ['\t'.join(x) for x in map]
        maps.append('')
        logf.write('%s %s: Writing map file\n' % (thisprog,timenow()))
        outmap.write('\n'.join(maps))
    else:
        head = f.next().strip()
        rslist = head.split()
    nrs = len(rslist) # number of markers
    elen = 2*nrs + 6 # expected # elements on each line
    logf.write('%s %s: found %d for nrs\n' % (thisprog,timenow(),nrs))
    eig = {}
    eig['founders'] = [array.array('c',[]) for x in xrange(nrs)] # marker rows, subject cols
    eig['offspring'] = [array.array('c',[]) for x in xrange(nrs)] # marker rows, subject cols
    adicts = [{} for x in xrange(nrs)] # count of alleles in a dict for each marker
    refallele = [None for x in xrange(nrs)] # list of first observed alleles
    nsubj = 0
    indiv = {'founders':[],'offspring':[]}
    for lnum,l in enumerate(f):
        ll = l.strip().split()
        if (lnum+1) % 200 == 0:
            logf.write('%s %s: Processing line %d\n' % (thisprog, timenow(),lnum+1))
        if len(ll) < elen: # ? short ?
            logf.write('%s %s: Line %d is %d long, expected %d\n' % (thisprog, timenow(),lnum,len(ll),elen))
        else:
            nsubj += 1
            sid = '%s_%s' % (ll[0],ll[1])
            isFounder = isOff = False
            status = labels[0] # founder unaff
            if ll[2] <> '0' and ll[3] <> '0': # has parent ids
               iclass = 'offspring'
               status = labels[2] # unaffected offspring
               if ll[5] == '2':
                    status = labels[3] # affected offspring
            else:
               iclass = 'founders'
               if ll[5] == '2':
                   status = labels[1] #change from unaff to aff founder label
            gender = 'M'
            if ll[4] == '2':
                gender = 'F'
            indiv[iclass].append('%s %s %s' % (sid,gender,status)) # for the ind file
            for snp in xrange(nrs):
                g1,g2 = ll[2*snp + 6],ll[2*snp + 7] # pair of genos
                if mdict.get(g1,None) or mdict.get(g2,None): # one or both missing
                    esnp = emissval # missing value
                else:
                    if not refallele[snp]:
                        refallele[snp] = g1 # first one we saw!
                    for g in (g1,g2):
                        n = adicts[snp].get(g,0)
                        n += 1
                        adicts[snp][g] = n
                    if g1 == g2: # hom
                        if g1 == refallele[snp]:
                            esnp = ehom2 # 2 copies of current reference allele
                        else:
                            esnp = ehom1 # no copies
                    else:
                        esnp = ehet # het - always has one copy of reference allele
                eig[iclass][snp].append(esnp) # append the eigenstrat geno code for this new subject
    for ek in eig.keys():
        lek = len(eig[ek])
        if len(eig[ek]) > 0:
            lek0 = len(eig[ek][0])
            s = 'for %s, have %d snp of len %d' % (ek,lek,lek0)
            print s
            logf.write(s)
            for x in range(lek):
                    if len(eig[ek][x]) <> lek0:
                        s = 'for row %d, len = %d, not %d' % (x, len(eig[ek][x]),lek0)
                        print s
                        logf.write(s)
    logf.write('%s %s: Now checking major allele assignment and fixing as needed\n' % (thisprog,timenow()))

    pfname = os.path.join(repository,basename,'.ped')
    mfname = os.path.join(repository,basename,'.map')
    try:
        pf = file(pfname,'r')
    except:
        pf = None
    try:
        mf = file(mfname,'r')
    except:
        mf = None
    if pf and mf: # have files will transpose from subject to snpwise
        
    
    
    return genos, snpAnno

def splitShareMatrix(chromlist=['22',],outFname='FramGainQC'):
    """challenge is to combine the various release categories
       and add their sample ids and columns on to grow the file?
    """
    outGaingeno = file('%s.geno' % outFname,'w')
    outGainmap = file('%s.snp' % outFname,'w')
    outGainped = file('%s.ped' % outFname,'w')
    sampleIds = []
    ped = getPed()
    p = ped.values()
    p.sort()
    pres = ['\t'.join(x) for x in p]
    outGainped.write('\n'.join(pres))
    outGainped.close()
    genos,snpAnno = getGenos()
    # now have all subjects in genos and markers in snp # append these                                
    a = snpAnno.values()
    print 'a=',a[:10]
    a.sort(lambda x,y: cmp((x[1],x[2]),(y[1],y[2]))) # in chrom offset order
    rslist = [x[0] for x in a] # rs in order

    print 'rslist=',rslist[:10]
    if chromn == 0: # first time only
        outGaingeno.write('markerID\t')
        outGaingeno.write('\t'.join(sampleIds))
        outGaingeno.write('\n')
        outGainmap.write('PREFERRED_ID\tCHR\tPOSITION\tQC_TYPE\tSEQUENCE\n')
    for rs in rslist:
        outGaingeno.write('%s\t' % rs)
        outGaingeno.write('\t'.join(genos[rs]))
        outGaingeno.write('\n')
        sl = snpAnno[rs]
        sl[2] = '%d' % sl[2]
        outGainmap.write('\t'.join(sl))
        outGainmap.write('\n')
    outGainmap.close()
    outGaingeno.close()          


if __name__ == "__main__":
        splitShareMatrix()   
