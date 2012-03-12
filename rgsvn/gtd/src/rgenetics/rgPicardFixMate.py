import os, sys
import shutil, tempfile
import optparse
import string
from subprocess import Popen

import pysam

from rgutils import timestamp
from samtools import Bam, Sam

progname = os.path.split(sys.argv[0])[1]

"""
fix mate pair information
Copyright ross lazarus oct 31 (happy halloween!)
All rights reserved
released for rgenetics under the LGPL
"""



class fixMate():
    def __init__(self, opts):
        if not os.path.exists(opts.jar): raise ValueError, 'jar file does not exist'
        self.jar = opts.jar
        self.tmp_dir = opts.tmp_dir
        self.log_suffix = opts.suffix_tmpl % (timestamp(), 'log')
        self.maxjheap = opts.maxjheap
        self.output_dir = os.path.dirname(opts.output)
        self.runs = []
        self.logs = []
        
    def run(self, input_file, output_file=None):
        assert input_file != None

        log_fd, log_filename = tempfile.mkstemp(dir=opts.tmp_dir, suffix=self.log_suffix)
        log_file = open(log_filename, 'w')

        # When no output file is defined, Fixmate tool will replace
        # the content of the input file.
        params = {'I=':input_file}
        if output_file: 
            params['O='] = output_file

        cmd_string = 'java -Xmx%s -jar %s %s' % (
            self.maxjheap, self.jar, ' '.join([''.join(it) for it in list(params.items())])
            )
        
        process = Popen(cmd_string, shell=True, stderr=log_file, stdout=log_file)
        
        return_value = process.wait()
        log_file.close()
        log_contents = open(log_filename).read()
        self.runs.append(cmd_string)
        self.logs.append(log_contents)
        
        shutil.move(log_filename, os.path.join(self.output_dir, os.path.basename(log_filename)))

        
if __name__ == '__main__':
    '''
    <command interpreter="python">
    rgPicardFixMate.py -i "$input_file" -o "$out_file" --tmp_dir "${__new_file_path__}"  
    --newformat "$newformat" -j "${GALAXY_DATA_INDEX_DIR}/shared/jars/FixMateInformation.jar"
   </command>
    '''
    op = optparse.OptionParser()
    op.add_option('-j','--jar')
    op.add_option('-i', '--input', default=None)
    op.add_option('--input-type', default=None)
    op.add_option('-o', '--output', default=None)
    op.add_option('--newformat', default='bam')
    op.add_option('-x', '--maxjheap', default='2g')
    op.add_option('--tmp_dir', default='/tmp')
    op.add_option('-n', '--title', default='')
    op.add_option('-v','--verbose', action='store_true')

    opts, args = op.parse_args()
    assert opts.input <> None
    assert os.path.isfile(opts.input) 
    assert os.path.isfile(opts.jar)

    if not os.path.exists(opts.tmp_dir): os.makedirs(opts.tmp_dir)

    suffix_tmpl = 'rgTempFixMate_%s.%s'
    opts.suffix_tmpl = suffix_tmpl
    
    # Picard tools always produces bam files. Naming accordingly...
    ifd, intermediate_filename = tempfile.mkstemp(suffix=suffix_tmpl % (timestamp(), 'bam'),
                                                  dir=opts.tmp_dir)

    fm = fixMate(opts)
    fm.run(os.path.abspath(opts.input), output_file=intermediate_filename)

    if opts.verbose:
        for run, log in zip(fm.runs, fm.logs):
            print 'Call to \n%s\n produced the following output:\n%s' % (run, log)

    intermediate = Bam(intermediate_filename)

    # Picard tool produced intermediate bam file. Depending on the
    # desired format, we either just move to final location or create
    # a sam version of it.
    if opts.newformat == 'sam':
        sam_file = intermediate.to_sam(os.path.abspath(opts.output))
        os.unlink(intermediate_filename)
    else:
        shutil.move(intermediate_filename, os.path.abspath(opts.output))
    
