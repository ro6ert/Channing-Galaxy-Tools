import os,string,sys,optparse,shutil
from subprocess import Popen
from rgutils import galhtmlprefix,galhtmlpostfix,galhtmlattr,timenow,fixPicardOutputs

progname = os.path.split(sys.argv[0])[1]

'''
picard wrapper
remember to set VALIDATION_STRINGENCY=LENIENT or you will go crazy like mike and I already did...

'''

if __name__ == '__main__':
    '''
    called as
    rgPicardGCBiasMetrics.py -i "$input_file" -d "$html_file.files_path" -o "$html_file"
     -w "$windowsize" -m "$mingenomefrac" -n "$out_prefix"
    -j ${GALAXY_DATA_INDEX_DIR}/shared/jars/CollectGCBiasMetrics.jar
#if $genomeSource.refGenomeSource == "history":
    -r "$genomeSource.ownFile"
#else:
    -r "$genomeSource.indices"
#end if

    called as
    java -jar /share/shared/galaxy/tool-data/shared/jars/CollectGcBiasMetrics.jar REFERENCE_SEQUENCE="hg18.fasta"
    MINIMUM_GENOME_FRACTION=0.00001 INPUT=test.bam OUTPUT=picardASMetrics.txt OUTPUT=test.txt CHART_OUTPUT=test.pdf
    WINDOW_SIZE=100 VALIDATION_STRINGENCY=LENIENT

    '''
    killme = string.punctuation + string.whitespace
    trantab = string.maketrans(killme,'_'*len(killme))
    op = optparse.OptionParser()
    op.add_option('-i', '--input', default=None)
    op.add_option('-o', '--html-output', default=None)
    op.add_option('-d', '--output-dir', default="/tmp/HsMetrics")
    op.add_option('-r', '--refseq', default='')
    op.add_option('-w', '--windowsize', default='100')
    op.add_option('-m', '--mingenofrac', default='0.00001')
    op.add_option('-j', '--jar', default='')
    op.add_option('-l', '--log', default=None)
    op.add_option('-n', '--namejob', default='GC Bias Metrics')
    op.add_option('-x', '--maxjheap', default='4g')
    op.add_option('--tmp_dir', default='/share/shared/tmp')
    opts, args = op.parse_args()
    try:
        os.makedirs(opts.output_dir)
    except:
        pass
    title = opts.namejob.translate(trantab)
    tempout = os.path.join(opts.output_dir,'rgPicardGCBiasMetrics.out')
    temppdf = os.path.join(opts.output_dir,'rgPicardGCBiasMetrics.pdf')
    temptab = os.path.join(opts.output_dir,'rgPicardGCBiasMetrics.xls')
    opts.log_file = opts.log or os.path.join(opts.output_dir, '%s.log' % title)
    # Create output folder and save our R script in there.
    standard_fd = open(opts.log_file, 'w')
    cl = ['java -Xmx',opts.maxjheap,' -jar ',opts.jar,' REFERENCE_SEQUENCE=',opts.refseq,' WINDOW_SIZE=',opts.windowsize,
    ' MINIMUM_GENOME_FRACTION=',opts.mingenofrac,' INPUT=',opts.input,' OUTPUT=',tempout,' TMP_DIR=',opts.tmp_dir,
    ' CHART_OUTPUT=',temppdf,' SUMMARY_OUTPUT=',temptab,' VALIDATION_STRINGENCY=LENIENT']
    process = Popen(''.join(cl), shell=True, stderr=standard_fd, stdout=standard_fd, cwd=opts.output_dir)
    return_value = process.wait()
    cl = 'mogrify -format jpg -resize x500 %s' % (temppdf) # make the jpg for fixPicardOutputs to find
    process = Popen(cl, shell=True, stderr=standard_fd, stdout=standard_fd, cwd=opts.output_dir)
    return_value = process.wait()
    standard_fd.close()
    fixPicardOutputs(tempout=temptab,output_dir=opts.output_dir,
      log_file=opts.log_file,html_output=opts.html_output,progname=progname,cl=cl)










