import os, sys
import string 
import shutil, optparse
from subprocess import Popen

from picard import bedToPicInterval, fixPicardOutputs

progname = os.path.split(sys.argv[0])[1]


'''
picard wrapper
remember to set VALIDATION_STRINGENCY=LENIENT or you will go crazy like mike and I already did...

'''

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
    op = optparse.OptionParser()
    op.add_option('-i', '--input', default=None)
    op.add_option('-o', '--html-output', default=None)
    op.add_option('-d', '--output-dir', default="/tmp/HsMetrics")
    op.add_option('-b', '--bait', default='')
    op.add_option('-t', '--target', default='')
    op.add_option('-j', '--jar', default='')
    op.add_option('-l', '--log', default=None)
    op.add_option('-n', '--namejob', default='HsMetrics')
    op.add_option('-x', '--maxjheap', default='2g')
    op.add_option('--tmp_dir', default='/export/tmp')
    opts, args = op.parse_args()
    try:
        os.makedirs(opts.output_dir)
    except:
        pass

    title = opts.namejob.translate(trantab)

    # Making sure that html_out has a name.
    if not opts.html_output: 
        opts.html_output = os.path.abspath(os.path.join(opts.output_dir, '%s.html' % title))

    # Convert bait and target files to bed format
    bait_file = os.path.join(opts.output_dir,'rgPicardHsMetrics.bait')
    target_file = os.path.join(opts.output_dir,'rgPicardHsMetrics.target')
    bedToPicInterval(infile=os.path.abspath(opts.bait), outfile=bait_file)
    bedToPicInterval(infile=os.path.abspath(opts.target), outfile=target_file)

    temp_out_file = os.path.join(opts.output_dir,'rgPicardHsMetrics.out')

    opts.log_file = opts.log or os.path.join(opts.output_dir, '%s.log' % title)
    sfd = open(opts.log_file, 'w')

    cmd = 'java -Xmx%s -jar %s BI=%s TI=%s I=%s O=%s VALIDATION_STRINGENCY=LENIENT TMP_DIR=%s' % (
        opts.maxjheap, os.path.abspath(opts.jar), 
        bait_file, target_file, os.path.abspath(opts.input), 
        temp_out_file, opts.tmp_dir)

    process = Popen(cmd, shell=True, stderr=sfd, stdout=sfd, cwd=opts.output_dir)
    return_value = process.wait()
    sfd.close()
    fixPicardOutputs(tempout=temp_out_file, output_dir=opts.output_dir, log_file=opts.log_file,
                     html_output=opts.html_output, progname=progname, cl=cmd)



    
 
    
    
