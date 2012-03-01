# read an eigenstrat eigenvalue file
# calculate euclidean distance and find closest pairs
# ross for rgGRR
# this code now avoids the horrors of R distance matrics

import math,sys

def calcMeanSD(useme):
    """
    A numerically stable algorithm is given below. It also computes the mean.
    This algorithm is due to Knuth,[1] who cites Welford.[2]
    n = 0
    mean = 0
    M2 = 0

    foreach x in data:
      n = n + 1
      delta = x - mean
      mean = mean + delta/n
      M2 = M2 + delta*(x - mean)      // This expression uses the new value of mean
    end for

    variance_n = M2/n
    variance = M2/(n - 1)
    """
    mean = 0.0
    M2 = 0.0
    sd = 0.0
    n = len(useme)
    if n > 1:
        for i,x in enumerate(useme):
            delta = x - mean
            mean = mean + delta/(i+1) # knuth uses n+=1 at start
            M2 = M2 + delta*(x - mean)      # This expression uses the new value of mean
        
        variance = M2/(n-1) # assume is sample so lose 1 DOF 
        sd = pow(variance,0.5)
    return mean,sd

def usingR(fname = 'c:/ross/rgGRR/mkkdist.xls'):
    """ deprecated in favour of going straight to the eigenstrat evec file

    test using R to calculate distance matrix and write it out
    popname='MKK'
    #mkk = read.table('c:/ross/rgGRR/MKK_founders_pca.xls.evec',head=F,sep='\t')
    eig = read.delim(paste('c:/ross/rgGRR/',popname,'_founders_pca.xls.evec',sep=''),head=F,sep = "", comment.char='~',skip=1,row.names=1)
    eigd = data.matrix(mkk[,1:10])
    d = dist(eigd)
    pdf(file=paste(popname,'dendro.pdf',sep=''),10,8)
    plot(hclust(d),main=paste(popname,"Pairwise Eigenvector Distance Dendrogram"),cex=0.5)
    grid(10,100)
    grid(10,5,col="maroon")
    abline(h=0, col='red')
    dev_.off()
    write.table(data.frame(as.matrix(d)),file=paste('c:/ross/rgGRR/',popname,'dist.xls',sep=''),sep='\t')    
    """

    f = file(fname,'r')
    colnames = None
    res = [] # list of (distance, id1,id2) tuples
    for row in f:
        if colnames == None:
            colnames = row.split() # note no rowname header!
        else:
            lrow = row.split()
            id = lrow.pop(0) # now aligned with colnames
            lrow = [(float(x),colnames[i]) for i,x in enumerate(lrow)]
            lrow.sort() # note that the diag will always be first
            dist,id2 = lrow[1] # should be the smallest nondiag
            res.append('%f %s %s' % (dist,id,id2))
    res.sort()
    print '\n'.join(res[:20])

def Euclidean(v1=[],v2=[]):
    """ euclidean distance
    e = sqrt(pow(p1-q1,2) + pow(p2-q2,2)...)
    """
    len1 = len(v1)
    assert len1 > 0
    assert len1 == len(v2)
    dist = math.sqrt(sum([math.pow(v1[i] - v2[i],2) for i in range(len1)]))
    return dist



def DistFromEvec(fname='c:/ross/rgGRR/MKK_founders_pca.xls.evec',maxTokeep=10):
    """ read an eigenstrat output eval file
    skip header, estimate all pairwise euclidean distances and
    save top 1% or so sorted into distance order
    The top pairs are likely to be highly related - check with rgGRR output eg
    """
    f = file(fname,'r').readlines()
    eigenvals = f.pop(0).split()[1:] # drop eigenvals row and '#eigvals:'
    dat = [x.split() for x in f] # matrix - first col is sample id
    statuses = [] # keep eigenstrat individual status - we use foundernonaff eg
    for i,row in enumerate(dat):
        for j in range(1,len(row)-1): # ignore status column
            row[j] = float(row[j])
        statuses.append(row.pop(-1)) # keep status
        dat[i] = row # drop status
    n = len(dat) # number of samples
    alldist = []
    distrib = []
    for i in range(n-1):
        id1 = dat[i][0]
        v1 = dat[i][1:] # ignore sample id
        distances = []
        for j in range(i+1,n):
            assert i <> j # for testing - no harm leaving in I guess - this seems pretty quick
            id2 = dat[j][0]
            v2 = dat[j][1:]
            dist = Euclidean(v1,v2)
            distances.append((dist,id2,id1))
            distrib.append(dist)
        distances.sort() # keep top n 
	alldist += distances[:maxTokeep] # only keep a subset
    alldist.sort() # so have closest at start
    nearest = ['%f\t%s\t%s' % (x[0],x[1],x[2]) for x in alldist] # ready for printing
    return nearest, distrib



if __name__ == "__main__":
   
    fname = sys.argv[1] # 'c:/ross/rgGRR/MKK_founders_pca.xls.evec'
    nearest,distrib = DistFromEvec(fname=fname)
    n = len(nearest)
    nd = len(distrib)
    top5pc = int(nd/20) # take 5pc of top pairs
    mean,sd = calcMeanSD(distrib)
    print '%d distances, mean = %g, sd = %g' % (nd,mean,sd)
    f = file('%s.distances' % fname,'w')
    s = ['%g' % x for x in distrib]
    f.write('\n'.join(s))
    f.write('\n')
    f.close()
    f = file('%s.toppairs' % fname,'w')
    f.write('\n'.join(nearest))
    f.write('\n')
    f.close()
    ntoprint = n/100 # top 1%
    s = ['%d\t%s' % (i,x) for i,x in enumerate(nearest[:ntoprint])]
    print '\nRanked minimum distances for pairs from %s' % fname
    print 'Rank\tDistance\tID1\tID2'
    print '\n'.join(s),'\n'


        
    
