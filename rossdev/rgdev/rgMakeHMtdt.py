# emergency!
# no useable tdt files?
# read the hapmap files and make some
# in all cases, make CEU offspring affected

def makeoffaff(f=None):
    """ return list of rows fiddled """
    fl = f.readlines()
    fl = [x.strip().split() for x in fl]
    for i,row in enumerate(fl):
        pa,ma = row[2:4]
        if pa <> '0' or ma <> '0':
            fl[i][5] = '2'
    return fl

def maketdt(orace='CEU',urace=['CHB','JPT'],chrom='22'):
    f = file('hm550k_chr%s_%s.ped' % (chrom,orace),'r')
    aff = makeoffaff(f)
    for u in urace:
        ul = file('hm550k_chr22_%s.ped' % u,'r').readlines()
        ul = [x.strip().split() for x in ul]
	aff += ul
    outf = 'hm%s_%s_chr%s_tdt.ped' % (orace,''.join(urace),chrom)
    f = file(outf,'w')
    saff = [' '.join(x) for x in aff]
    f.write('\n'.join(saff))
    f.write('\n')
    f.close()
    print 'wrote %d lines to %s' % (len(aff),outf)
    rl = ['%d:%d' % (i,len(row)) for i,row in enumerate(aff)]
    print ','.join(rl)

if __name__ == "__main__":
    maketdt()
    maketdt(urace=['YRI'])

