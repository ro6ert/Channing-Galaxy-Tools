# january 2008 for rgenetics/rgalaxy
# copyright 2008 ross lazarus
# note the soe files have multiple events to ignored

import sys,os,glob,time

fdatapath="/nfs/framinghamShare_nov2007/pheno" # where the unrar'd things be

def getPhes(listglob=None,minsub=1000):
    """try pfts
    """
    res = {}
    lsacc = []
    ltacc = []
    lphv = []
    allvars = {}
    badrows = []
    repevents = []

    def isFE(cname=''):
        """ don't want main questionnaire things at present
        """
        noint = 0
        if cname[:2] == 'FE' or cname[:2] == 'MF': # ignore most main exam things
            try:
                n = int(cname[2:])
            except:
                n = 0
            if n:
                return 1
        else:
            #print 'not isFE(%s)' % cname
            return 0
            
    def addFPheno(thisglob="*.pft*",fdatapath="/nfs/framinghamShare_nov2007/pheno"):
        """
        amalgamate all data by shareid for each study type
        Be sure to pass an unambiguous fglob
        files have typically:
        
        [rerla@hg pheno]$ head -20 *pft*_1s.*txt | more
        ==> phs000040.pht000102.v1.p1.c1.pft3_1s.GRU.txt <==
        # Study accession: phs000040.v1.p1
        # Table accession: pht000102.v1.p1.c1
        # Consent group: General Research Use
        # Citation instructions: The study accession (phs000040.v1.p1) is used to cite the study and its data tables and documents. The data in this 
        file should
         be cited using the accession pht000102.v1.p1.c1.
        # To cite columns of data within this file, please use the variable (phv#) accessions below:
        #
        # 1) the table name and the variable (phv#) accessions below; or
        # 2) you may cite a variable as phv#.v1.p1.c1.

        ##phv00022567   phv00022568     phv00022569     phv00022570     phv00022571     phv00022572     phv00022573     phv00022574     phv00022575     
        phv00022
        576     phv00022577     phv00022578     phv00022579
        IDTYPE  na_1_3  fv1_1_3 fv6_1_3 fvc_1_3 pf_1_3  ft_1_3  rat_1_3 mmf_1_3 ff1_1_3 ff2_1_3 ff3_1_3 shareid
        3       3       3.259363413     4.529647827     4.750863075     8.383769035     14.20197868     0.686   2.340326786     4.004788399     
        3.514954329
        0.963076949     5
        3       3       3.904792309     4.417870522     4.466477871     8.011305809     7.794268131     0.874   5.51417017      7.987002373     
        6.561185837
        2.777911425     11

        we need the data and the variable names,id's and table/study accessions
        """
        ffglob = os.path.join(fdatapath,thisglob)
        flist = glob.glob(ffglob)
        print '# Note: flist=',flist
        sid = 'shareid' # need to find this to merge all data
        for fname in flist:
            sa = 'Study'
            ta = 'Table'
            saccession = None
            taccession = None
            f = file(fname,'r')
            hfound = 0
            head = None
            phv = None
            cshareid = None
            distinguish = fname.split('.')[-3]
            nshort = 0
            colswewant = []
            for n,line in enumerate(f):
                if hfound == 0:
                    if line.find(sa) == 2:
                        saccession = line.split(': ')[1] # get study id
                        lsacc.append(saccession)
                    elif line.find(ta) == 2:
                        taccession = line.split(': ')[1] # get table
                        ltacc.append(taccession)
                    elif line[:2] == '##':
                        phv = line[2:].split() # list of variable ids
                        lphv.append(phv)
                        hfound = 1
                elif hfound == 1:
                    hfound = 2
                    head = line.strip().split('\t') # real variable names
                    try:
                        cshareid = head.index(sid)
                    except:
                        print '### problem in file %s - cannot find shareid column in header %s' % (fname,head)
                        sys.exit(1)
                    colswewant.append(cshareid)
                    minl = len(head)
                    if fglob == '*.ex*':
                        for c,cname in enumerate(head):
                           if cname == 'FE68':
                               cname = 'heightIn'
                               head[c] = cname
                               colswewant.append(c)
                           elif cname == 'FE67':
                               cname = 'weightLb'
                               head[c] = cname
                               colswewant.append(c)
                    else:
                        for c,cname in enumerate(head):
                            colswewant.append(c) # keep all
                    #head = ['%s_%s' % (x.upper(),distinguish) for x in head]
                    #head[cshareid] = sid # eesh
                    for col in colswewant:
                        cname = head[col]
                        if not allvars.get(cname,None):
                            allvars[cname] = 0
                else:
                    row = line.strip().split('\t')
                    if len(row) >= minl:
                        id = row[cshareid] # get shareid we hope
                        sofar = res.get(id,{}) # new dict if not seen yet
                        for c in colswewant:
                            cname = head[c]
                            if sofar.get(cname,None):
                                if row[c] <> sofar[cname]:
                                    print '## Error - c %d row %d shareid %s file %s has %s for %s - was %s' % (c,
                                            n,id,fname,row[c],cname,sofar[cname])
                            elif row[c] > '':
                                sofar[cname] = row[c]
                                allvars[cname] += 1
                            else:
                                sofar[cname] = 'NA'
                        res[id] = sofar
                    else:
                        nshort += 1
            if nshort > 0:
                s = 'Errors: file %s has %d cases short for headlen %d' % (fname,nshort,minl)
                badrows.append(s)


    for thisGlob in listglob:
        addFPheno(thisglob=thisGlob)
        try:
            firstrec = res.values()[0]
        except:
            firstrec = {}
        cols = len(firstrec.keys())
        print 'thisGlob=%s. res has %d cols' % (thisGlob,cols)
    if len(badrows) > 0:
        print '\n'.join(badrows)
        print ''
    if len(repevents) > 0:
        print '\n'.join(repevents)
        print ''
    s = res.keys() # shareids
    s = [int(x) for x in s]
    s.sort()
    s = ['%d' % x for x in s]
    vk = [x for x in allvars.keys() if allvars[x] > minsub and allvars[x] <> 'shareid']
    vk.sort()
    vk.insert(0,'shareid')
    ivk = [x for x in allvars.keys() if allvars[x] <= minsub]
    ivk.sort()
    if len(ivk) > 0:
        ires = ['%s:%d ' % (x,allvars[x]) for x in ivk]
        print '###! ignoring phenotypes with fewer than %d values' % minsub
        print '###! %s' % ','.join(ires)
    outfname = 'framSharePheno.xls'
    outf = file(outfname,'w')
    outf.write('%s\n' % '\t'.join(vk))
    print 'outputing',vk
    for id in s:
        row = [res[id].get(v,'NA') for v in vk]
        row = '\t'.join(row)
        outf.write('%s\n' % row)
    outf.close()
    print 'wrote %d rows and %d cols to %s' % (len(s),len(vk),outfname)
    
if __name__ == "__main__":
    if len(sys.argv) > 1:
        listglob = sys.argv[1].split(',')
    else:
        listglob = ['*.pft*','*.ex*','*.vitd*','*.tnf*','*adi*','*abdomct*','*bsc*']
    if len(sys.argv) > 2:
        minsub = int(sys.argv[2])
    else:
        minsub = 1
    print 'using listglob=%s and minimum subjects per pheno = %d' % (listglob,minsub)
    getPhes(listglob=listglob,minsub=minsub)

