# copyright 2010 ross lazarus
# released under the LGPL
#

import os
import sys
import datetime
import string
import tempfile
from subprocess import Popen

from mako.template import Template

from rgutils import getFileString, timenow

class PicardTool(object):
    WRAPPER_BASENAME = 'UnnamedPicardTool'
    HTML_REPORT_TEMPLATE = 'templates/picard_reports.mako'
    
    def run(self, cl=None, redir=True):
        assert cl <> None
        
        fd, templog = tempfile.mkstemp(dir=self.opts.outdir, suffix=self.__class__.WRAPPER_BASENAME) 
        tlf = open(templog,'w')
        if redir:
           process = Popen(' '.join(cl), shell=True, stderr=tlf, stdout=tlf, cwd=self.opts.outdir)
        else:
           process = Popen(' '.join(cl), shell=True, stderr=tlf, cwd=self.opts.outdir)
        return_value = process.wait()
        tlf.close()
        tlogs = ''.join(open(templog,'r').readlines()).strip()
        return_msg = '## executed %s; returning STATUS %d.' % (' '.join(cl), return_value)
        log_msg = tlogs or 'Nothing appeared on stderr/stdout'

        os.unlink(templog) # always

        self.tlog.write('\n'.join([return_msg, 'Log records:', log_msg]))

        return return_value
        
    def write_html_report(self):
        """
        write the report as html
        """
        
        out_folder = self.opts.outdir
        my_template = Template(
            filename=os.path.join(self.tool_folder, self.__class__.HTML_REPORT_TEMPLATE), 
            strict_undefined=True)
        
        if os.path.exists(os.path.abspath(out_folder)):
            files = [os.path.join(out_folder, x) for x in os.listdir(out_folder) 
                     if not x.startswith('.')]
            files.sort(key=lambda f: os.path.getmtime(f))
            file_info = [(os.path.split(f)[-1], getFileString(os.path.split(f)[-1], self.opts.outdir))
                         for f in files]
        else:
            file_info = []
              
        template_parameters = {
            'program_name': str(self.program_name), 
            'timestamp' : str(timenow()),
            'file_info': file_info,
            'log_data' : open(self.tlogname).read().replace('\n', '\n<BR />'),
            'command_string': str(self.cl)
            }
            
        f = open(self.opts.htmlout, 'w')
        f.write(my_template.render(**template_parameters))
        f.close()

        


    def cleanup(self):
        if self.tidy and hasattr(self, 'delme'):
            for fname in self.delme:
                try:
                    os.unlink(fname)
                except:
                    pass
        s = open(self.tlogname,'r').readlines()
        print >> sys.stdout, s # for info
        os.unlink(self.tlogname)


 


