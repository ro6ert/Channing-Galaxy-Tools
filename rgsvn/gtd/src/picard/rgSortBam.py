import os,string,sys,optparse,shutil,tempfile
from subprocess import Popen
from rgutils import galhtmlprefix,galhtmlpostfix,galhtmlattr,timenow,keynat
progname = os.path.split(sys.argv[0])[1]
tidy = False # delete intermediate temp files or keep when testing
verbose = True

"""
From http://www.broadinstitute.org/gsa/wiki/index.php/Input_files_for_the_GATK
- Do we always need to do this after bwa has generated an alignment?
Eeesh..

Fixing BAM files with alternative sortings

The GATK requires that the BAM file be sorted in the same order as the
reference. Unfortunately, many BAM files have headers that are sorted
in some other order -- lexicographical order is a common alternative.
To resort the BAM file so that it matches your reference, you will
have to strip out the BAM file's reference, replacing it with the
reference from your header. To do this:
Convert the bam -> sam, stripping away the header
Reimport the bam, using your reference's fai as the ref_index. This
step results in an 'unsorted' bam with a correctly sorted header.
Sort the bam file using 'samtools sort'.
Index the file
For example:
samtools view <your file>.bam > -o <new file>.sam
samtools import <your reference>.fai <your file>.sam <your
file>.sorted_header.bam
samtools sort <your file>.sorted_header.bam <your file>.sorted
samtools index <your file>.sorted.bam

galaxy users will need to option to replace headers using a built in or history fasta index
so much stuffing around needed below to preserve any existing metadata
and to index a history fasta if needed

"""
def sortBAM(opts=None):
    """
    work broken out here
    """

    def getSamExtraHeaders(tempsam=None):
        """
        """
        ehead = []
        f = open(tempsam,'r')
        for row in f:
            if not row.startswith('@'): # end of headers
                break
            if not row.startswith('@SQ'):
                ehead.append(row)
        f.close()
        return ehead

    def run(cl=None,tlog=None):
        process = Popen(' '.join(cl), shell=True, stderr=tlog, stdout=tlog)
        rval = process.wait()
        if verbose:
            print 'executing %s returned %d' % (' '.join(cl),rval)

    def makeFaidx(opts=None,delme=[]):
        # must make a temp softlink to index as samtools faidx is so dumb it will try to write [inputpath].fai 
        fd,tempfasta = tempfile.mkstemp(dir=opts.tmp_dir,suffix='rgfakeFasta.fasta')
        os.unlink(tempfasta) # we just want the name
        delme.append(tempfasta)
        newsamfai = '%s.fai' % tempfasta
        delme.append(newsamfai)
        cl = ['ln -s',opts.fasta,tempfasta]
        process = Popen(' '.join(cl), shell=True, stderr=tlog, stdout=tlog)
        rval = process.wait()
        cl = ['samtools faidx',opts.fasta]
        process = Popen(' '.join(cl), shell=True, stderr=tlog, stdout=tlog)
        rval = process.wait()
        if verbose:
            print 'executing %s returned %d' % (' '.join(cl),rval)
        return newsamfai
        
    def viewSam(infile=None):
        """
        """
        fd,tempsam = tempfile.mkstemp(dir=opts.tmp_dir,suffix='rgSortBamTemp.sam')
        cl = ['samtools view -h -o',tempsam,infile]
        if verbose:
            print 'convert %s to sam' % infile
        process = Popen(' '.join(cl), shell=True, stderr=tlog, stdout=tlog)
        rval = process.wait()
        if verbose:
            print 'executing %s returned %d' % (' '.join(cl),rval)
        return tempsam
    
    def runCalmd(infile=None):
        """
        samtools calmd assuming a bam - convenient to do it with sortedbam before adding rest of header metadata
        """
        fd,tempbam = tempfile.mkstemp(dir=opts.tmp_dir,suffix='rgCalmdTemp.bam')
        if opts.fai: # punt that we have a built  
            refseq = os.path.splitext(opts.fai)[0] # ugh
        else:
            refseq = opts.fasta
        cl = ['samtools calmd -b',infile,refseq,'>',tempbam]
        run(cl,tlog)
        cl = ['mv',tempbam,infile] # replace what came in destructive but
        run(cl,tlog)
        
    def runCleanSam(infile=None,isSam=True):
        """
        optional trim of reads beyond refseq
        """
        fd,tempsam = tempfile.mkstemp(dir=opts.tmp_dir,suffix='rgCleanSam.sam')
        if isSam:        
            cl = ['java -Xmx2g -jar',opts.cleansam,'INPUT=%s' % infile,'OUTPUT=%s',tempsam]
            run(cl,tlog)
            cl = ['mv',tempsam,infile] # replace what came in destructive but
            run(cl,tlog)
        else: # given bam we have to do some work
            sam = viewSam(infile=infile)
            cl = ['java -Xmx2g -jar',opts.cleansam,'INPUT=%s' % sam,'OUTPUT=%s',tempsam]
            run(cl,tlog)
            cl = ['samtools view -h -b -S -o',infile,tempsam] # replace what came in destructive but
            run(cl,tlog)             
                    
    delme = []
    killme = string.punctuation + string.whitespace
    trantab = string.maketrans(killme,'_'*len(killme))
    title = opts.title.translate(trantab)
    fd,tlogname = tempfile.mkstemp(dir=opts.tmp_dir,suffix='rgSortBam.log')
    tlog = open(tlogname,'w') 
    delme.append(tlogname)
    fd,sortedbam = tempfile.mkstemp(dir=opts.tmp_dir,suffix='rgSortBamSorted.bam') 
    delme.append(sortedbam)
    fd,newbam = tempfile.mkstemp(dir=opts.tmp_dir,suffix='rgSortBamnewbam.bam')
    delme.append(newbam) # if we need bam, we can copy this to galaxy output
    if opts.informat == 'bam': # need text so use samtools view
        tempsam = viewSam(infile=opts.input)
    else:
        tempsam = opts.input # is already sam
    ehead = getSamExtraHeaders(tempsam=tempsam)
    if len(ehead) == 0:
        print >> sys.stdout, '## ehead is empty on %s' % tempsam
    if opts.fasta > '': # supplied genome fasta needs indexing
        newsamfai = makeFaidx(opts=opts,delme=delme)    
    else: # use supplied fai
        newsamfai = opts.fai
    cl = ['samtools import',newsamfai,tempsam,newbam] # stitch the adjusted headers
    run(cl,tlog)
    cl = ['samtools sort',newbam,sortedbam]
    run(cl,tlog)
    sortedbam = '%s.bam' % sortedbam # samtools insists on adding .bam
    if opts.calmd: # use the -b flag and remember to redirect back into sortedbam
        runCalmd(infile=sortedbam)
    # sortedbam is go
    if len(ehead) > 0: # eesh - need to view and add the additional non SQ metadata headers back again
        if verbose:
            print'## WARNING putting back non @SQ metadata (%s) from %s' % (ehead,opts.input)
        fd,tempout = tempfile.mkstemp(dir=opts.tmp_dir,suffix='rgfixheadout.sam')
        delme.append(tempout)
        tempsam = viewSam(infile=sortedbam)
        delme.append(tempsam)
        f = open(tempout,'w') # keep repaired file here
        f.write(''.join(ehead))
        for row in open(tempsam,'r').readlines():
            f.write(row)
        f.close()
        # fixed stuff is in tempout
        if opts.cleansam > '':
            runCleanSam(infile=tempsam,isSam=True) # since we have a sam, do it now
        if opts.newformat == 'bam':
            cl = ['samtools view -h -b -S -o',opts.output,tempout]
            run(cl,tlog)
        else:
            cl = ['cp',tempout,opts.output]
    else: 
        if opts.cleansam > '':
            runCleanSam(infile=sortedbam,isSam=False) # sigh. TODO: some extra conversions could be avoided with more code
        if opts.newformat == 'bam':    
            cl = ['cp',sortedbam,opts.output]
        else:
            tempsam = viewSam(infile=sortedbam)
            delme.append(tempsam)
            cl = ['cp',tempsam,opts.output]
        run(cl,tlog)
    if tidy:
        for f in delme:
            os.unlink(f) # cleanup

    
if __name__ == '__main__':
    '''
   <command interpreter="python">
   rgSortBam.py -i "$input_file" --informat "$input_file.ext" -o "$out_file" -n "$out_prefix" 
   --newformat "$new_format" --tmp_dir "${__new_file_path__}"
#if $runCleanSam:
--cleanSam
#end if
#if $runCalmd:
--calmd
#end if
#if $newHead.refGenomeSource=="indexed":
   --fai "$newHead.indexsrc"
#else
   --fasta "$newHead.indexsrc"
#end if
   </command>
    '''
    op = optparse.OptionParser()
    op.add_option('-i', '--input', default=None)
    op.add_option('-o', '--output', default=None)
    op.add_option('-n', '--title', default="SortSamBam")
    op.add_option('--newformat', default='bam')
    op.add_option('--informat', default='bam')
    op.add_option('-x', '--maxjheap', default='2g')
    op.add_option('--tmp_dir', default='/tmp')
    op.add_option('--fai', default='')
    op.add_option('--fasta', default='')
    op.add_option('--calmd', default=False)
    op.add_option('--cleansam', default='') # will have the cleansam jar path
    opts, args = op.parse_args()
    assert opts.input <> None
    assert os.path.isfile(opts.input) 
    assert not (opts.fasta > '' and opts.fai > '')
    try:
        os.makedirs(opts.tmp_dir)
    except:
        pass
    sortBAM(opts=opts)
    
    
    
