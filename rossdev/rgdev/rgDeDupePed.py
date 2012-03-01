# not so quick fix for a sharpe problem
# grr does not like duplicate samples
# or subjects with only one parent
# so this script cleans em out
# and optionally takes a random subset of markers
# 2500 is enough usually
# june 2008 ross lazarus

import sys,random,os,copy

progname = os.path.split(sys.argv[0])[-1]

def main(infname='test.ped',n=5000,fulltrios=1):
    """can't run abecasis' grr with dupes or missing parents
    simple filter to clean those records out - takes first of dupes only

    """
    keepme = None # signal not to take a random subset of snps
    f = file(infname,'r')
    if n:
        outfname = 'nodupe_%d_%s' % (n,infname)
    else:
        outfname = 'nodupe_%s' % infname
    o = file(outfname,'w')
    allids = {}
    mas = {}
    pas = {}
    idks = {}
    # first pass - record mother and father
    for i,row in enumerate(f):
        if i == 0:
            nmarkers = (len(row.strip().split()) - 6)/2
        if i % 1000 == 0 and i <> 0:
            print 'first pass at line',i
        fid,iid = row.strip().split()[:2]
        k = '%s_%s' % (fid,iid)
        allids[k] = k

    # try to make a map file to match
    mfname = '%s.map' % os.path.splitext(infname)[0]
    human = ['%d' % x for x in range(1,23)] # plink uses 23..26 or so
    # so we do too...
    humanautos = dict(zip(human,human)) # for quick lookup
    try:
        maprows = file(mfname,'r').readlines()
        autosomal = [i for i,x in enumerate(maprows) if humanautos.get(x[0],None)]
        nmap = len(autosomal) # number of markers
        nsubset = n
    except:
        maprows = None
        nsubset = n
        nmap = n
        autosomal = range(nmap)
        print '### cannot find map %s matching %s - no subset map written' % (mfname,infname)
    if n <> 0: # must subset
        keepme = autosomal # make a list of marker offsets
        random.shuffle(keepme)
        keepme = keepme[:nsubset] # fast way to get a random subset
        keepme.sort() # keep in original mapfile order
        print 'using %d long keepme - %s' % (len(keepme),keepme[:20])
    else:
        nmap = nmarkers
        keepme = None
    if maprows: # rewrite
        newmapname = 'nodupe_%d_%s' % (nsubset,mfname)
        if keepme:
            newmap = [maprows[x] for x in keepme]
            s = 'random'
        else:
            newmap = maprows # write entire file
            s = '(all)'
        newmapf = file(newmapname,'w')
        newmapf.write(''.join(newmap))
        newmapf.close()
        print '## wrote %d %s markers to new map file %s' % (len(keepme),s,newmapname)



    if nmap <> nmarkers:
                print '### strange - map %s has %d rows, but first row of %s has %d markers - WTF?' % (mfname,
                    nmap,infname,nmarkers)
    f.seek(0) # rewind for second pass
    nrecs = nout = mmiss = fmiss = dupes = 0
    for i,row in enumerate(f):
        nrecs += 1
        if i % 1000 == 0 and i <> 0:
            print 'second pass at line',i
        fid,iid,faid,maid = row.strip().split()[:4]
        k = '%s_%s' % (fid,iid)
        if idks.get(k,None):
            dupes += 1
            print '## skipping duplicate for fid %s iid %s' % (fid,iid)
            continue
        idks[k] = k
        if faid <> '0' or maid <> '0':
            fk = '%s_%s' % (fid,faid)
            mk = '%s_%s' % (fid,maid)
            ik = '%s_%s' % (fid,iid)
            if fulltrios and not allids.get(fk,None):
                print '# Father %s of %s not found - skipping' % (fk,ik)
                fmiss += 1
                continue # reject - missing parent
            if fulltrios and not allids.get(mk,None):
                print '# Mother %s of %s not found - skipping' % (mk,ik)
                mmiss += 1
                continue # reject - missing parent
        if keepme: # must subset the row using snp in keepme
            lrow = row.strip().split()
            newrow = lrow[:(6 + 2*nsubset)] # snarf ped and a list of the right size
            for n,col in enumerate(keepme): # slow
                outoffset = 2*n + 6 # first allele to write in newrow
                inoffset = 2*col + 6 # where to look for first allele in lrow
                newrow[outoffset] = lrow[inoffset]
                newrow[outoffset+1] = lrow[inoffset+1]
            o.write(' '.join(newrow))
            o.write('\n')
        else: # no subset
            o.write(row)
        nout += 1
    o.close()

    print '# removed %d subjects total of %d, %d missing fathers, %d missing mothers, %d dups' % (nrecs-nout,
                                                        nrecs,fmiss,mmiss,dupes)


def usage(n=1):
    """print instructions
    """
    print 'Remove duplicate subjects, optionally take a random subset of an lped file (write a map) and optionally'
    print 'remove any offspring from incomplete trios (missing parents) to keep GRR happy'
    print 'python %s input_pedfile_path [n markers to subset randomly] [f to remove incomplete trio offspring]' % progname
    sys.exit(n)

if __name__ == "__main__":
    if len(sys.argv) > 1:
        infname = sys.argv[1]
    else:
        print '## error - need an input ped file path as the CL parameter please'
        usage(1)
    if len(sys.argv) > 2:
        n = sys.argv[2]
        try:
            n = int(n)
            if n > 0:
                print '## will subset to %d markers' % n
        except:
            print '## error - if you want a random subset of markers, second CL param must be an integer count - eg 5000'
            usage(2)
    else:
        n = 0 # no need for subset
    if len(sys.argv) > 3:
        fulltrios = sys.argv[3].lower() == 'f'
    else:
        fulltrios = 0
    main(infname=infname,n=n, fulltrios = fulltrios)

>>>>>>> .r1364
