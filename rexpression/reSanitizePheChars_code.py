#build list of available data
import string, os, sys, glob, shutil, logging
import galaxy.util
from galaxy import datatypes, config, jobs 

gal_Log = logging.getLogger(__name__)

repository = "/usr/local/galaxy/data/rg/library"


builds = []
#Read build names and keys from galaxy.util
for dbkey, build_name in galaxy.util.dbnames:
    builds.append((build_name,dbkey,False))

def get_available_builds(defval='hg18'):
    for i,x in enumerate(builds):
        if x[1] == defval:
           x = list(x)
           x[2] = True
           builds[i] = tuple(x)
    return builds


# for caco
# called as  plinkCaCo.py $i $name $test $outformat $out_file1 $logf $map
from galaxy import datatypes
import time,string

def timenow():
    """return current time as a string
    """
    return time.strftime('%d/%m/%Y %H:%M:%S', time.localtime(time.time()))

def exec_after_process(app, inp_data, out_data, param_dict, tool, stdout, stderr):
    """Sets the name of the data"""
    job_name = param_dict.get( 'title', 'Sanitize' )
    killme = string.punctuation + string.whitespace
    trantab = string.maketrans(killme,'_'*len(killme))
    title = job_name.translate(trantab)
    n,indat = inp_data.items()[0]
    new_ext = indat.ext
    gal_Log.debug('## execafter - ext = %s' % new_ext)
    data = out_data.items()[0][1]
    data = app.datatypes_registry.change_datatype(data, new_ext)
    data.name = job_name
    app.model.context.flush()
    
