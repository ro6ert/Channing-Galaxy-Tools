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
    run this Picard tool as
    java -jar /share/shared/galaxy/tool-data/shared/jars/CollectAlignmentSummaryMetrics.jar REFERENCE_SEQUENCE="hg18.fasta" 
    ASSUME_SORTED=true ADAPTER_SEQUENCE='' IS_BISULFITE_SEQUENCED=false MAX_INSERT_SIZE=100000 INPUT=test.bam 
    OUTPUT=picardASMetrics.txt VALIDATION_STRINGENCY=LENIENT

    called as
        rgPicardASMetrics.py -i "$input_file" -d "$html_file.files_path" -o "$html_file"
    -s "$sorted" -b "$bisulphite" -a "$adaptors" -m $maxinsert -n "$out_prefix"
    -j ${GALAXY_DATA_INDEX_DIR}/shared/jars/CollectAlignmentSummaryMetrics.jar
#if $genomeSource.refGenomeSource == "history":
    -r "$genomeSource.ownFile"
#else:
    -r "$genomeSource.indices"
#end if

    '''
    killme = string.punctuation + string.whitespace
    trantab = string.maketrans(killme,'_'*len(killme))
    op = optparse.OptionParser()
    op.add_option('-i', '--input', default=None)
    op.add_option('-o', '--html-output', default=None)
    op.add_option('-d', '--output-dir', default="/tmp/HsMetrics")
    op.add_option('-r', '--refseq', default='')
    op.add_option('-s', '--assume-sorted', default='true')
    op.add_option('-a', '--adapters', action='append', type="string")
    op.add_option('-b', '--bisulphite', default='false')
    op.add_option('-m', '--maxinsert', default='100000')
    op.add_option('-j', '--jar', default='')
    op.add_option('-l', '--log', default=None)
    op.add_option('-n', '--namejob', default='Alignment Summary Metrics')
    op.add_option('-x', '--maxjheap', default='2g')
    op.add_option('--tmp_dir', default='/tmp')
    opts, args = op.parse_args()
    try:
        os.makedirs(opts.output_dir)
    except:
        pass
    title = opts.namejob.translate(trantab)
    tempout = os.path.join(opts.output_dir,'rgPicardASMetrics.out')
    opts.log_file = opts.log or os.path.join(opts.output_dir, '%s.log' % title)
    # Create output folder and save our R script in there.
    standard_fd = open(opts.log_file, 'w')
    cl = ['java -Xmx',opts.maxjheap,' -jar ',opts.jar,' REFERENCE_SEQUENCE=',opts.refseq,' ASSUME_SORTED=',
        opts.assume_sorted,''.join([' ADAPTER_SEQUENCE=%s' % x for x in opts.adapters]), 
        ' IS_BISULFITE_SEQUENCED=',opts.bisulphite, ' MAX_INSERT_SIZE=',opts.maxinsert,
        ' INPUT=',opts.input,' OUTPUT=',tempout,' VALIDATION_STRINGENCY=LENIENT',' TMP_DIR=',opts.tmp_dir]
    process = Popen(''.join(cl), shell=True, stderr=standard_fd, stdout=standard_fd, cwd=opts.output_dir)
    return_value = process.wait()
    standard_fd.close()
    
    fixPicardOutputs(tempout=tempout,output_dir=opts.output_dir,log_file=opts.log_file,html_output=opts.html_output,
       progname=progname,cl=cl)




        


    

    
    
