import os,string,sys,optparse,shutil
from subprocess import Popen
from rgutils import galhtmlprefix,galhtmlpostfix,galhtmlattr,timenow,bedToPicInterval,fixPicardOutputs

progname = os.path.split(sys.argv[0])[1]


'''
picard wrapper
remember to set VALIDATION_STRINGENCY=LENIENT or you will go crazy like mike and I already did...

'''
# note base.css has a table.colored but it sucks

usage_msg = """
python rgPicardHsMetrics.py -i $input_file -d $html_file.files_path -o $html_file
    -b $bait_bed -t $target_bed -n "$out_prefix" --tmp_dir "${__new_file_path__}"
    -j ${GALAXY_DATA_INDEX_DIR}/shared/jars/CalculateHsMetrics.jar  
"""



if __name__ == '__main__':
    '''
    run as
    java -jar /share/shared/galaxy/tool-data/shared/jars/CalculateHsMetrics.jar 
    BAIT_INTERVALS=test.pic TARGET_INTERVALS=test.pic INPUT=test.bam OUTPUT=picardHsMetrics.txt

    called as <command interpreter="python">
    rgPicardHsMetrics.py -i $input_file -d $html_file.files_path -o $html_file
    -b $bait_bed -t $target_bed -n "$out_prefix" --tmp_dir "${__new_file_path__}"
    -j ${GALAXY_DATA_INDEX_DIR}/shared/jars/CalculateHsMetrics.jar 
  </command>
    '''
    killme = string.punctuation + string.whitespace
    trantab = string.maketrans(killme,'_'*len(killme))
    op = optparse.OptionParser(usage=usage_msg)
    op.add_option('-i', '--input')
    op.add_option('-o', '--html-output')
    op.add_option('-d', '--output-dir', default="/tmp/HsMetrics")
    op.add_option('-b', '--bait')
    op.add_option('-t', '--target')
    op.add_option('-j', '--jar', default='')
    op.add_option('-l', '--log', default='')
    op.add_option('-n', '--namejob', default='HsMetrics')
    op.add_option('-x', '--maxjheap', default='2g')
    op.add_option('--tmp_dir', default='/tmp')
    opts, args = op.parse_args()

    has_inputs = opts.input and opts.bait and opts.target
    has_outputs = opts.html_output and opts.output_dir
    
    if not (has_inputs and has_outputs):
        print usage_msg
        sys.exit(-1)


    try:
        os.makedirs(opts.output_dir)
    except:
        pass
    title = opts.namejob.translate(trantab)
    baitf = os.path.join(opts.output_dir,'rgPicardHsMetrics.bait')
    targetf = os.path.join(opts.output_dir,'rgPicardHsMetrics.target')

    bedToPicInterval(opts.bait, baitf)

    if opts.target == opts.bait: # same file sometimes
        targetf = baitf
    else:
        bedToPicInterval(opts.target, targetf)   
    tempout = os.path.join(opts.output_dir,'rgPicardHsMetrics.out')
    opts.log_file = opts.log or os.path.join(opts.output_dir, '%s.log' % title)


    sfd = open(opts.log_file, 'w')

    cl = ['java -Xmx',opts.maxjheap,' -jar ',opts.jar,
          ' BAIT_INTERVALS=',baitf,' TARGET_INTERVALS=',
          targetf,' INPUT=',os.path.abspath(opts.input),
          ' OUTPUT=',tempout,
          ' VALIDATION_STRINGENCY=LENIENT', ' TMP_DIR=',opts.tmp_dir
          ]
    process = Popen(''.join(cl), shell=True, stderr=sfd, stdout=sfd, cwd=opts.output_dir)
    return_value = process.wait()
    sfd.close()

    fixPicardOutputs(tempout=tempout, output_dir=opts.output_dir, 
                     log_file=opts.log_file,
                     html_output=os.path.join(opts.output_dir, opts.html_output),
                     progname=progname,
                     cl=cl)



    

    
    
