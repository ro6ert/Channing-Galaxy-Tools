"""parse an fbat log
qqplot
"""

import sys,os, math

def parsefbatlog(fname=''):
    """look for lines of the format
    read in: 9 markers from 900 pedigrees (901 nuclear families,2757 persons)
    medelian error: locus rs1860184_731, pedigree 136 [1400,318]
    and results like
    A total of 3750286 mendelian errors have been found        
genotypes of families with mendelian error have been reset to 0
>> load share_pulmpheno.phe
6 quantitative traits have been successfully read          
9340 persons have been phenotyped
warning: 1191 persons are not in any pedigrees            
>> trait ppfev1adj
affection ppfev1adj** ppfvcadj ppfefadj ppratioadj ppfefratadj
>> fbat
trait ppfev1adj; offset 0.000; model additive; test bi-allelic; minsize 10; min_freq 0.000; p 1.000; maxcmh 1000

Marker        Allele   afreq   fam#     S-E(S)      Var(S)       Z           P
------------------------------------------------------------------------------
rs10900309    2        0.184    865     14.574     509.676   0.646    0.518558
rs10900309    4        0.816    865    -14.574     509.676  -0.646    0.518558
rs12039954    2        0.082    433     -5.929     204.286  -0.415    0.678284
    """
    logf = '%s.log' % fname
    mendelerrors = {}
    markers = {}
    res = []
    try:
        f = open(logf,'r')
    except:
        print '**No fbat log file %s found - Mendelian errors not removed' % logf
        return mendelerrors
    for l in f:
        if l.find('locus') <> -1: # might be of interest
            marker = ''
            pedigree = ''
            l = l.replace(',',' ')
            ll = l.strip().split()
            for i in range(len(ll)):
                if ll[i].lower() == 'locus': # warning - could change with fbat upgrades?
                    marker = ll[i+1]
                elif ll[i].lower() == 'pedigree':
                    pedigree = ll[i+1]
            if marker > '':
                if not markers.has_key(marker):
                    markers[marker] = []
                if not mendelerrors.has_key(pedigree):
                    mendelerrors[pedigree] = [marker,]
                else:
                    mendelerrors[pedigree].append(marker)
                if pedigree not in markers[marker]:
                    markers[marker].append(pedigree)
        elif l[:6].lower == 'marker':
            head = l.strip().split()
        elif l[:2] == 'rs': # results line
            ll = l.strip().split()
            row = [ll[n] for n in [0,7]]
            res.append(ll)
            #marker,allele,afreq,fam,ses,vars,z,p = ll
    return mendelerrors,res,markers

def runfbat(fbatpedf=None, smendelfile=None, mmendelfile=None, ggoutfile=None,
            phefile=None, phelist=[], centered=[], empirical=[]):
    """Swap fbat log files to get the mendel errors out, then get the results for each trait out separately for
    keeping, and amalgamated into a genome graphs file. Use a bunch of temp files and operate fbat in batch
    mode by feeding a list of commands from stdin as part of the Popen run...
    """
    fbate = "/usr/local/bin/fbat"
    clist = [] # these will be fed via stdin to fbat!
    clist.append('log %s' % (mendelfile)) # write all mendel errors here
    s = "load %s" % fbatpedf
    clist.append(s)
    traitlogs = []
    for i, phe in enumerate(phelist): # add trait and fbat flags to list
        (f,afilename) = tempfile.mkstemp(prefix=phe)
        traitlogs.append(afilename) # to use later to parse logs
        clist.append('trait %s' % phe)
        clist.append('log %s' % (afilename)) # switch log files - we'll parse them later
        task = ['fbat',]
        if len(c) > i:
            c = centered[i]
            e = empirical[i]
            if e:
                task.append('-e')
            if c:
                task.append('-o')
        s = ' '.join(task)
        clist.append(s)
    x = subprocess.Popen(fbate,shell=True,stdin=clist) # ignore stdout
    retval = x.wait()
    me,res,markers = parsefbatlog(mendelfile)
    mlist = []
    for m in markers.keys():
        mlist.append((len(markers[m]),m)) # replace with count of pedigrees, marker tuple for sorting
    mlist.sort()
    mlist.reverse()
    plist = []
    for p in mendelerrors.keys():
        plist.append((len(mendelerrors[p]),p))
    plist.sort()
    plist.reverse()
    res = []
    res.append('Marker\tN')
    outf = open(mmendelfile,'w')
    for (n,m) in mlist:
        res.append('%d\t%s' % (n,m))
    outf.close()
    outf = open(smendelfile,'w')    
    res.append('N\tPedigree')
    for (n,p) in plist:
        res.append('%d\t%s' % (n,p))
    res.append('')
    outf.write(string.join(res,'\n'))
    outf.close()
    fout = file(ggoutfile,'w')
    fout.write('marker\ttrait\t\tp\tlog10p\n'
    for n,afilename in enumerate(traitlogs):
        trait = phelist[n]
        me,res = parsefbatlog(afilename)
        # run fbat, then parse the output file
        for row in res:
            try:
                lp = -math.log10(row[1])
            except:
                lp = 1.0
            row.insert(1,trait)
            row.append('%e' % lp)
            fout.write('\t'.join(row))
            fout.write('\n')
    fout.close()        


def parse():
    """ 
    <command interpreter="python2.4">
        fbat.py $pedf $phef "$title" "$phecols" $out_file1 $logfile $mendelfile
    </command>
     """
    prog = os.path.split(sys.argv[0])[-1]
    print >> sys.stdout,'## Rgenetics http://rgenetics.org Galaxy Tools - cl=%s \n' % (sys.argv)
    pedfname = sys.argv[1]
    phefname = sys.argv[2]
    jobname = sys.argv[3]
    phecols = sys.argv[4].split()
    outfile = sys.argv[5]
    logfile = sys.argv[6]
    mendelfile = sys.argv[7]
    #centered=sys.argv[8].split()
    #empirical=sys.argv[9].split()
    runfbat(fbatpedf=pedfname, mendelfile=mendelfile, outfile=outfile,
            phefile=phefname, phelist=phecols, centered=[], empirical=[])    



if __name__ == "__main__":
    parse()