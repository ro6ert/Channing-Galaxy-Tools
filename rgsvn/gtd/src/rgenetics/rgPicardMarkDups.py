""" 
Mark duplicates with Picard
Copyright ross lazarus nov 17 2010
All rights reserved released for rgenetics under the LGPL 
"""

import os, sys
import string
import optparse
from subprocess import Popen 

from picard import PicardTool

progname = os.path.split(sys.argv[0])[1]

class markDups(PicardTool):

    WRAPPER_BASENAME = 'rgPicardMarkDups'

    def __init__(self, opts=None, cl=[], tidy=False):
        self.opts = opts
        self.tidy = tidy
        self.cl = ' '.join(cl) # ready for the htmlfile output
        self.delme = []
        killme = string.punctuation + string.whitespace
        trantab = string.maketrans(killme,'_'*len(killme))
        self.title = self.opts.title.translate(trantab)
        self.tlogname = os.path.join(self.opts.outdir, '%s.txt' % self.__class__.WRAPPER_BASENAME)
        self.tlog = open(self.tlogname, 'w')
        self.metricstxt = 'rgPicardMarkDupsMetrics.txt'

        program = os.path.abspath(cl[0])
        self.tool_folder = os.path.dirname(program)
        self.program_name = os.path.basename(program)
    
    
    def run_tool(self):
        """
        """

        cl = ['java -Xmx%s' % self.opts.maxjheap,'-jar', self.opts.jar, 
              'I=', self.opts.input, 'O=',self.opts.outbam,
              'M=', os.path.abspath(os.path.join(self.opts.outdir, self.metricstxt))
              ]
        if self.opts.remdups == 'true':
            cl.append('REMOVE_DUPLICATES=true')
        if self.opts.assumesorted == 'false':
            cl.append('ASSUME_SORTED=false')
        if self.opts.readregex > '':
            cl.append('READ_NAME_REGEX="%s"' % self.opts.readregex)
        if float(self.opts.optdupedist) > 0:
            cl.append('OPTICAL_DUPLICATE_PIXEL_DISTANCE=%s' % self.opts.optdupedist)
        
        run = self.run(cl)
        self.tlog.close()
           
 
    
if __name__ == '__main__':
    '''
      <command interpreter="python">
   rgPicardMarkDups.py -i "$input_file" -n "$out_prefix" --tmp_dir "${__new_file_path__}" 
#if $remdupes:
  --remdups "true"
#end if
    --assumesorted "$assumeSorted" --readregex "$readRegex" --optdupedist "$optDupeDist"
   -j "${GALAXY_DATA_INDEX_DIR}/shared/jars/MarkDuplicates.jar" -d "$html_file.files_path" -t "$html_file"
  </command>


    '''
    op = optparse.OptionParser()
    op.add_option('-i', '--input', default=None)
    op.add_option('-o', '--outbam', default=None)
    op.add_option('-n', '--title', default="markDupes")
    op.add_option('--remdups', default='false') 
    op.add_option('--assumesorted', default='true') 
    op.add_option('--readregex', default="[a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*")
    op.add_option('--optdupedist', default="100")
    op.add_option('-t', '--htmlout', default="")
    op.add_option('-d', '--outdir', default="")
    op.add_option('-x', '--maxjheap', default='2g')
    op.add_option('--tmp_dir', default='/tmp')
    op.add_option('-j','--jar',default='')
    opts, args = op.parse_args()
    assert opts.input <> None
    assert opts.outbam <> None
    assert os.path.isfile(opts.input)
    assert os.path.isfile(opts.jar)
    if not os.path.exists(opts.tmp_dir): os.makedirs(opts.tmp_dir)
    if not os.path.exists(opts.outdir): os.makedirs(opts.outdir)

    tool = markDups(opts=opts,cl=sys.argv)
    tool.run_tool()
    tool.write_html_report()
