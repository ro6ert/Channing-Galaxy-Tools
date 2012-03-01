#!/usr/bin/env python2.4

"""
Creates an html file to be used for viewing rgenetics files

   <command interpreter="python2.4">rgenetics_import.py $file_type_dir $base_name $output </command>


"""

import sys, os, glob, shutil,subprocess

galhtmlprefix = """<?xml version="1.0" encoding="utf-8" ?>
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
<head>
<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
<meta name="generator" content="Galaxy %s tool output - see http://g2.trac.bx.psu.edu/" />
<title></title>
<link rel="stylesheet" href="/static/style/base.css" type="text/css" />
</head>
<body>
<div class="document">
"""

def doProjectImport():
    """ import into an html composite data type for Rgenetics
        Important metadata manipulations in an exec_after_process hook:
    """
    progname = sys.argv[0]
    file_type_dirs = sys.argv[1].strip().split(',')
    base_names = sys.argv[2].strip().split(',')
    output = sys.argv[3]
    out = open(output,'w')
    out.write(galhtmlprefix % progname)
    out.write('<h2>Rgenetics data imported from %s</h2><hr/>\n<ul>' % ','.join(basenames))
    for n,file_type_dir in enumerate(file_type_dirs):
        basename = base_names[n]
        flist = glob.glob("%s.*" % os.path.join(file_type_dir, basename))
        for i, data in enumerate( flist ):
            out.write('<li><a href="%s">%s</a></li>\n' % (os.path.split(data)[-1],os.path.split(data)[-1]))
    out.write("</ul></div></body></html>")
    out.close()



def doImport():
    """ 
    import into one of the new html composite data types for Rgenetics
    Important metadata manipulations in an exec_after_process hook:        
    now called as   
  <command interpreter="python">rgenetics_import.py "$file_type" "$base_name" "$output" "$output.files_path" "${GALAXY_DATA_INDEX_DIR}"</command>
    </command>

    """
    progname = sys.argv[0]
    file_type = sys.argv[1]
    base_names = sys.argv[2]
    bnlist = base_names.split(' ')
    base_name = bnlist[0]
    outf = sys.argv[3]
    ofdir = sys.argv[4]
    file_lib = sys.argv[5]
    file_type_dir = os.path.join(file_lib,file_type)
    srcglob = "%s.*" % os.path.join(file_type_dir, base_name)
    out = open(outf,'w')
    out.write(galhtmlprefix % progname)
    try:
        os.makedirs(ofdir)
        out.write('## made %s\n' % ofdir)
    except:
        pass 
    out.write('<h2>Rgenetics data imported from %s</h2>\n<hr><ul>' % base_name)
    flist = glob.glob(srcglob)
    print '#!!@@ for %s, get %s' % (srcglob,flist)
    out.write('for %s, get %s' % (srcglob,flist))
    for i, fname in enumerate( flist ):
        print 'writing fname %s' % fname
        fbase = os.path.split(fname)[-1]
        out.write('<li><a href="%s">%s</a></li>\n' % (fbase,fbase))
        newhome = os.path.join(ofdir,fbase)
        shutil.copy(fname,newhome)
    out.write("</ul></div></body></html>")
    out.close()


if __name__ == "__main__": 
   doImport()

