import os, sys
import shutil, tempfile
import optparse
import string
from subprocess import Popen

import pysam

from rgutils import fixPicardOutputs

progname = os.path.split(sys.argv[0])[1]
verbose = True
tidy = False

'''
picard wrapper
for validatesam tool

'''

def samToBam(infile=None,outdir=None):
    """
    use samtools view to convert sam to bam
    """
    fd,tempbam = tempfile.mkstemp(dir=outdir,suffix='rgSortBamTemp.bam')
    cl = ['samtools view -h -b -S -o',tempbam,infile]
    self.run(cl)
    return tempbam

def cleanSam(insam=None, newsam=None, picardErrors=[],outformat=None,sortme=True):
    """
    stub - need to do the work of removing all the error sequences
    pysam is cool
    infile = pysam.Samfile( "-", "r" )
    outfile = pysam.Samfile( "-", "w", template = infile )
    for s in infile: outfile.write(s)
    
    errors from ValidateSameFile.jar look like
    WARNING: Record 32, Read name SRR006041.1202260, NM tag (nucleotide differences) is missing
    ERROR: Record 33, Read name SRR006041.1042721, Empty sequence dictionary.
    ERROR: Record 33, Read name SRR006041.1042721, RG ID on SAMRecord not found in header: SRR006041

    """
    removeNames = [x.split(',')[1].replace(' Read name ','') for x in picardErrors if len(x.split(',')) > 2]
    remDict = dict(zip(removeNames,range(len(removeNames))))
    infile = pysam.Samfile(insam,'r')
    print 'found %d identifiable error sequences in picardErrors out of %d' % (len(removeNames),len(picardErrors))
    if len(removeNames) > 0:
        outfile = pysam.Samfile(newsam,'w')
        i = 0
        for s in infile:
            if not remDict.get(s.qname,None): # not to be cleaned
                outfile.write(s)
            else:
                i += 1
        print 'Wrote %s after removing %d rows from %s' % (newsam,i,insam)
        outfile.close()
        infile.close()
    else: # we really want a nullop or a simple pointer copy
        infile.close()
        if newsam:
            shutil.copy(insam,newsam)
    
    

def picVal(opts=None):
    """
    called with sam so no need to convert
    """
    assert opts <> None
    killme = string.punctuation + string.whitespace
    trantab = string.maketrans(killme,'_'*len(killme))
    title = opts.title.translate(trantab)
    tempout = os.path.join(opts.output_dir,'rgPicardValidate.out') 
    temptab = os.path.join(opts.output_dir,'rgPicardValidate.xls')
    opts.log_file = opts.log or os.path.join(opts.output_dir, 'rgPicardValidate_%s.log' % title)
    # Create output folder and save our R script in there.
    stf = open(opts.log_file, 'w')
    sortedfile=None
    if verbose:
        print '# opts.ignore',opts.ignore,' opts.sortme=',opts.sortme
    if opts.sortme:
        fd,sortedfile = tempfile.mkstemp(suffix='rgcleansam.sorted.bam')
        if opts.datatype == 'sam': # need to work with a bam 
            tempbam = samToBam(opts.input,opts.outdir)
            pysam.sort(tempbam,sortedfile)
        else: # is already bam
            pysam.sort(opts.input,sortedfile)
    cl = ['java -Xmx',opts.maxjheap,' -jar ',opts.jar,' O=',tempout,' TMP_DIR=',opts.tmp_dir]
    if verbose:
        print '# cl so far',cl
    if opts.sortme:
    	cl.append(' I=%s' % sortedfile)
    else:
    	cl.append(' I=%s' % opts.input)
    if opts.maxoutput == '0':
        opts.maxoutput = '65535'
    cl.append(' MAX_OUTPUT=%s' % opts.maxoutput)
    if opts.ignore[0] <> 'None': # picard error values to ignore
        cl += [' IGNORE=%s' % x for x in opts.ignore if x <> 'None']
    if opts.bisulphite.lower() <> 'false':
        cl.append(' IS_BISULFITE_SEQUENCED=true')
    if opts.refseq <> '':
        cl += [' R=%s' % opts.refseq,]
    s1 = ' '.join(['"%s"' % x for x in cl])
    s = '## rgPicardValidate.py about to Popen:\n%s\n' % s1    
    stf.write(s)
    if verbose:
        print s    
    pefilename = os.path.join(opts.output_dir, 'rgPicardValidate_%s.errors' % title)
    picerrors = open(pefilename,'w')
    process = Popen(''.join(cl), shell=True, stderr=picerrors, stdout=picerrors, cwd=opts.output_dir)
    return_value = process.wait()
    picerrors.close()
    pe = open(pefilename,'r').readlines()
    stf.write('## got %d rows - first few =%s\n' % (len(pe), '\n'.join(pe[:5])))
    if opts.dryrun <> 'dryrun': # want to run cleansam
        if opts.dryrun == 'sam':
            outformat = 'sam'
            newsam = opts.sam
        elif opts.dryrun == 'bam':
            outformat = 'bam'            
            newsam = opts.bam
        cleanSam(insam=opts.input, newsam=newsam, picardErrors=pe,outformat=outformat,sortme=opts.sortme)
    stf.close()
    fixPicardOutputs(tempout=tempout,output_dir=opts.output_dir,
      log_file=opts.log_file,html_output=opts.html_output,progname=progname,cl=cl,transpose=False)
    if opts.sortme:
        os.unlink(sortedfile)
        if opts.datatype == 'sam': # was converted 
            os.unlink(tempbam) # temporary

if __name__ == '__main__':
    '''
    called as
     <command interpreter="python">
    rgPicardValidate.py -i "$input_file" --datatype "$input_file.ext" -d "$html_file.files_path" -o "$html_file"
    -t "$out_prefix" -e "$ignore" -b "$bisulphite" -m "$maxerrors" -y "$new_format"
    -j ${GALAXY_DATA_INDEX_DIR}/shared/jars/ValidateSamFile.jar
#if $genomeSource.refGenomeSource == "history":
    -r "$genomeSource.ownFile"
#elif $genomeSource.refGenomeSource=="indexed":
    -r "$genomeSource.indices"
#end if
#if $new_format=='sam':
 --sam "$out_file"
#elif $new_format=='bam':
 --bam "$out_file"
#end if
  </command>

    '''
    op = optparse.OptionParser()
    op.add_option('-i', '--input', default=None)
    op.add_option('-o', '--html-output', default=None)
    op.add_option('-d', '--output-dir', default="/tmp/validateSam") 
    op.add_option('-r', '--refseq', default='')
    op.add_option('-e', '--ignore', action='append', type="string")
    op.add_option('-y', '--dryrun', default='dryrun')
    op.add_option('--sam', default='None')
    op.add_option('--bam', default='None')
    op.add_option('-m', '--maxoutput', default='0')
    op.add_option('-b', '--bisulphite', default='false')
    op.add_option('-j', '--jar', default='')
    op.add_option('-l', '--log', default='rgValidateSam.log')
    op.add_option('-t', '--title', default='Validate Sam/Bam')
    op.add_option('-x', '--maxjheap', default='8g')
    op.add_option('-s', '--sortme', default=False)
    op.add_option('--tmp_dir', default='/tmp')
    op.add_option('--datatype', default='bam')
    opts, args = op.parse_args()
    try:
        os.makedirs(opts.output_dir)
    except:
        pass
    picVal(opts = opts)
       
