"""
## added genetic model tab delim file
## may 2008 rml
## released under the LGPL for the rgenetics/rgalaxy project
## changed to use a generator to create single output lines
## one at a time - the framingham data easily outstrips the size of
## a python string if I use the usual lazy everything at once approach
## autoconvert - run everything that's not there
## TODO make more selective and configurable
## class oct 2007 ross lazarus
## changed to use arrays rather than lists
## for snp vectors to save some space and maybe some time
## this is slow compared to convertf but doesn't require
## a whole slew of configuration parameters in a text file :)
## linkage format ped file to eigenstrat
## ross lazarus me fecit august 14 2007
## convertf seems to want a 7 column pedigree
## eigenstrat format is transposed with 1 line per snp and
## number of ref alleles for each marker for each subject
## subjects are in indiv file as sampleid gender label
## markers are in a linkage pedigree map file format
## how to avoid holding entire ped file in ram?
## read a line,
## need an allele count
## start by assuming first seen allele is major, then reverse if wrong
"""

import sys,os,array,time,string
from optparse import OptionParser

thisprog = os.path.basename(sys.argv[0])
additive = 'additive'
dominant = 'dominant'
recessive = 'recessive'
zadditive = 'zadditive'

progname = os.path.split(sys.argv[0])[1]

def timenow():
    """return current time as a string
    """
    return time.strftime('%d/%m/%Y %H:%M:%S', time.localtime(time.time()))


def readMap(fprefix="example"):
    """
    """
    map = []
    f = file('%s.map' % fprefix,'r')
    for l in f:
        ll = l.strip().split()
        if len(ll) >= 4:
            map.append((ll[:4]))
    return map

def majAllele(adict={'A':20,'C':100}):
    """return commonest allele
    """
    ak = adict.keys()
    if len(ak) == 1: # mono
        return ak[0] # is major allele
    m = None
    maxn = -9
    for k,n in adict.iteritems():
        if n > maxn:
            maxn = n
            m = k
    return m


class pedGenModelConverter:
    """ additive = # rare alleles 0,1,2 or -1,0,1?
    dom = 0,1,1 recessive = 0,0,1
    read a ped file and map file.

    Separate out founders and offspring since might only want one or the other

    Genotypes go into arrays to save ram, code them by guessing the
    major allele, and then fixing if we guessed wrong once we've figured out the allele
    frequencies.
    These are written out recoded into a genetic model for glm or other downstream
    uses, as a kind of excel compatible tab delimited
    recoded ped file - header row with column names and data
    """

    def __init__(self,pedfname = None,title='title',outfname='rgMakeGenModel.phe',
                 xtradir="/foo/bar", linkage=True, logf = None, model=additive,
                 labels=['FounderNonAff','FounderAff','OffspringNonAff','OffspringAff']):
        """ convert a linkage ped file into an excel style tab delimited file
        with allele pairs converted into codes representing the required genetic model
        """
        ehet = '1'
        ehom1 = '0'
        ehom2 = '2'
        missval = '.'
        missing = ['N','0','.','-']
        if model == additive: # we have ehom2 = 2 major alleles
            modeldict = {ehom1:'2',ehet:'1',ehom2:'0',missval:missval}
        elif model == zadditive:
            modeldict = {ehom1:'1',ehet:'0',ehom2:'-1',missval:missval}
        elif model == dominant:
            modeldict = {ehom1:'1',ehet:'1',ehom2:'0',missval:missval}
        else:
            modeldict = {ehom1:'2',ehet:'0',ehom2:'0',missval:missval}
        # model code translation dicts

        swapdict = {ehom1:ehom2,ehom2:ehom1,ehet:ehet,missval:missval} # if initial ref allele was wrong, use these to swap
        mdict = dict(zip(missing,missing))
        pedf = '%s.ped' % pedfname
        f = file(pedf,'r')
        if linkage: # read map file
            map = readMap(pedfname)
            rslist = [x[1] for x in map] # get rs numbers
        else:
            head = f.next().strip()
            rslist = head.split()
        nrs = len(rslist) # number of markers
        elen = 2*nrs + 6 # expected # elements on each line
        logf.write('%s %s: found %d for nrs\n' % (thisprog,timenow(),nrs))
        gm = {}
        gm['founders'] = [] # array.array('c',[]) for x in xrange(nrs)] # marker rows, subject cols
        gm['offspring'] = [] # [array.array('c',[]) for x in xrange(nrs)] # marker rows, subject cols
        adicts = [{} for x in xrange(nrs)] # count of alleles in a dict for each marker
        refallele = [None for x in xrange(nrs)] # list of first observed alleles
        nsubj = 0
        indiv = {'founders':[],'offspring':[]}
        for lnum,l in enumerate(f):
            ll = l.strip().split()
            if (lnum+1) % 200 == 0:
                logf.write('%s %s: Processing line %d\n' % (thisprog, timenow(),lnum+1))
            if len(ll) < elen: # ? short ?
                logf.write('%s %s: Line %d is %d long, expected %d\n' % (thisprog, timenow(),
                                            lnum,len(ll),elen))
            else:
                garray = array.array('c',[missval,]*nrs)
                nsubj += 1
                sid = '%s_%s' % (ll[0],ll[1])
                if sid == '1_1': # eesh
                   sid = '%d_%d' % (nsubj,nsubj)
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
                ped = ll[:6]
                indiv[iclass].append(ped) # for row wise output
                for snp in xrange(nrs):
                    pos = 2*snp + 6 # first
                    g1,g2 = ll[pos],ll[pos+1] # pair of genos
                    if mdict.get(g1,None) or mdict.get(g2,None): # one or both missing
                        esnp = missval # missing value
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
                    garray[snp] = esnp
                gm[iclass].append(garray) # append genos for this new subject
        for ek in gm.keys():
            lek = len(gm[ek])
            if len(gm[ek]) > 0:
                lek0 = len(gm[ek][0])
                s = 'for %s, have %d subjects with %d markers' % (ek,lek,lek0)
                print s
                logf.write(s)
                for x in range(lek):
                        if len(gm[ek][x]) <> lek0:
                            s = 'for row %d, len = %d, not %d' % (x, len(gm[ek][x]),lek0)
                            print s
                            logf.write(s)
        logf.write('%s %s: Now checking major allele assignment and fixing as needed\n' % (thisprog,timenow()))
        for iclass in gm.keys():
            for subject in xrange(len(gm[iclass])): # for each subject
                for snp in xrange(nrs): # now check to see if reference = major allele
                    major = majAllele(adicts[snp])
                    if major <> refallele[snp]: # either None or we need to change all the codes
                        if major <> None:
                            gm[iclass][subject][snp] = swapdict[gm[iclass][subject][snp]]
        for iclass in gm.keys():
            for subject in xrange(len(gm[iclass])): # now convert to genetic model wanted
                gm[iclass][subject] = [modeldict[x] for x in gm[iclass][subject]] # translate to model
        self.gm = gm
        self.model = model
        self.indiv = indiv
        self.nrs = nrs
        self.nsubj = nsubj
        self.logf = logf
        self.basename = title
        self.outfname = outfname
        self.rslist = rslist
        pedhead = 'famid\tiid\tfid\tmid\tgender\taffection\t'
        rs = '\t'.join(rslist)
        self.outhead = '%s%s\n' % (pedhead,rs)

    def gmRowGenerator(self,outclasslist=['founders']):
        """
        not working transposed for this output
        """
        for iclass in outclasslist: # have to write fully row wise
            for subject in xrange(len(self.gm[iclass])):
                row = self.indiv[iclass][subject] + list(self.gm[iclass][subject])
                # pedigree and genetic model codes for each snp - not a native rgenetics
                # type yet
                s = '\t'.join(row)
                yield s # we're a generator so iterable. Neat.


    def writeOut(self,fo=True,oo=False,):
        """Write the converted file remember to append cases to each snp row
        make everything tab delim for generic uses
        """
        s = 'writeOut fo = %s,oo = %s' % (fo,oo)
        # print s
        self.logf.write(s)
        nfound = len(self.indiv['founders'])
        noff = len(self.indiv['offspring'])
        if fo and oo: # all
            outclass = ['founders','offspring']
        elif oo: # do offspring only
                outclass = ['offspring',]
        elif fo: #  founders only
                outclass = ['founders']
        outf = file(self.outfname,'w')
        outf.write(self.outhead) # header row with pedigree and marker names
        res = self.gmRowGenerator(outclasslist=outclass)
        for x in res:
            # yes, this is slower than outf.write('\n'.join(res)) but that fails with framingham
            # because the file is bigger than a python string allows!
            outf.write('%s\n' % x)
        outf.close()



def makeGenModel(title='rgMakeGenModel', outfname='rgMakeGenModel.phe',
                 logf=None, pedfile=None, model='additive',
                 founders=True, offspring=False):
    """

    """
    pedfname = os.path.splitext(pedfile)[0] # get root
    mappath = '%s.map' % pedfname
    pedpath = '%s.ped' % pedfname
    pedf = file(pedpath,'r')
    test1 = len(pedf.readline().strip().split())
    test2 = len(pedf.readline().strip().split())
    pedf.close()
    if test1 <> test2: # must be fbat
        linkage = False
    else:
        linkage = True
    s = '## rgMakeGenModel: title=%s, pedpath=%s,cl=%s\n' % (title,
            pedfname, sys.argv)
    print >> sys.stdout, s # so will appear as blurb for file
    logf.write(s)
    c = pedGenModelConverter(pedfname=pedfname,title=title, outfname=outfname,
            linkage=linkage,logf = logf, model=model)
    nfound = len(c.indiv['founders'])
    noff = len(c.indiv['offspring'])
    c.writeOut(fo = founders, oo = offspring)



u = """rgMakeGenModel -o ./foo.html -p /usr/local/galaxy/data/rg/library/lped/foo.ped -m a -x N """


if __name__ == "__main__":
    """
    Updated oct 4 - converters now need to clone the behaviour of creating a new
    rgenetics data file which is really an html file


    <command interpreter="python2.4">
    rgPedEigConv.py $i.extra_files_path/$i.metadata.base_name $i.metadata.base_name
    "$title" $outfile.files_path $outfile $founders $offspring "$uf" "$af" "$uo" "$ao"
    </command>

    Pass a path to lped and an optional force rebuild flag
    glob and make everything if force else only if missing

    """
    ss = '%s%s' % (string.punctuation,string.whitespace)
    ptran =  string.maketrans(ss,'_'*len(ss))
    parser = OptionParser(usage=u, version="%s 1.2" % progname)
    a = parser.add_option
    a("-t","--title",type="str",dest='title',default='rgMakeGenModel',
      help="Descriptive title for outputs")
    a("-m","--gmodel",type="str",dest='gmodel',default='a',
      help="Genetic model for Model files - a=additive, d=dominant, r=recessive")
    a("-p","--pedfile",type="str",dest='pedfile',default=None,
      help="Input linkage format pedigree file")
    a("-o","--outfile",type="str",dest='ofname',default='rgConvert.phe',
      help="Output path for genetic model phe file")
    a("-l","--logfile",type="str",dest='lfname',default='rgConvert.log',
      help="Log file path")
    a("-f","--founders",action="store_true",dest='founders',
      help="Include founders")
    a("-c","--offspring",action="store_true",dest='offspring',
      help="Include offspring")
    (options,args) = parser.parse_args()
    gmod = "additive"
    if options.gmodel.lower() == 'r':
        gmod = "recessive"
    elif options.gmodel.lower() == 'd':
        gmod = "dominant"
    logfile = file(options.lfname,'w')
    title = options.title.translate(ptran) # for outputs
    makeGenModel(title=title,outfname=options.ofname,logf=logfile,
        pedfile=options.pedfile,model=gmod,
        founders=options.founders,offspring=options.offspring)
    try:
        lf.close()
    except:
        pass



