import os,string,sys,optparse,shutil
import shutil
import pdb
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
    op.add_option('-d', '--output-dir', default="/tmp/AsMetrics")
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

    # We must be able to create the output folder.
    output_folder = os.path.abspath(opts.output_dir)
    if not os.path.exists(output_folder): os.makedirs(output_folder)


    # Opts.adapters needs to be one string. If it's None, it must be coerced to empty str.
    opts.adapters = opts.adapters or ''

    title = opts.namejob.translate(trantab)
    tempout = os.path.join(output_folder,'rgPicardASMetrics.out')
    opts.log_file = opts.log or os.path.join(opts.output_dir, '%s.log' % title)


    # Picard only accepts reference files that have a 'fasta' extension name
    if not opts.refseq.endswith('fasta'):
        new_filename = '%s.fasta' % opts.refseq
        destination = os.path.join(output_folder, new_filename)
        shutil.copy(os.path.abspath(opts.refseq), destination)
        opts.refseq = destination
                    

    standard_fd = open(opts.log_file, 'w')
    cl = ['java -Xmx',opts.maxjheap,' -jar ',opts.jar,' REFERENCE_SEQUENCE=', os.path.abspath(opts.refseq),' ASSUME_SORTED=',
        opts.assume_sorted,''.join([' ADAPTER_SEQUENCE=%s' % x for x in opts.adapters]), 
        ' IS_BISULFITE_SEQUENCED=',opts.bisulphite, ' MAX_INSERT_SIZE=',opts.maxinsert,
        ' INPUT=',os.path.abspath(opts.input),' OUTPUT=',tempout,' VALIDATION_STRINGENCY=LENIENT',' TMP_DIR=',opts.tmp_dir]
    process = Popen(''.join(cl), shell=True, stderr=standard_fd, stdout=standard_fd, cwd=opts.output_dir)
    return_value = process.wait()


    standard_fd.close()

    fixPicardOutputs(tempout=tempout, output_dir=opts.output_dir, 
                     log_file=opts.log_file,
                     html_output=os.path.join(opts.output_dir, opts.html_output),
                     progname=progname,cl=cl)




        


    

    
    
