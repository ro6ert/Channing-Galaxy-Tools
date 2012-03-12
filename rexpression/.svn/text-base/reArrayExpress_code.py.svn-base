#build list of available data
import string, os, sys, glob, shutil
import galaxy.util
from galaxy import datatypes, config, jobs 


repository = "/usr/local/galaxy/data/rg/library"


galhtmlprefix = """<?xml version="1.0" encoding="utf-8" ?>
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
<head>
<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
<meta name="generator" content="Galaxy %s tool output - see http://galaxyproject.org/" />
<title></title>
<link rel="stylesheet" href="/static/style/base.css" type="text/css" />
</head>
<body>
<div class="document">
<ol>
"""



galhtmlpostfix = """</ol></div></body></html>"""

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


#post processing, set build for data and add additional data to history
def exec_after_process(app, inp_data, out_data, param_dict, tool, stdout, stderr):
    """tricky. We wrote     
    outlist = ['### makenewgalaxy\t%s\t%s\t%s' % (outfiles[i],outtypes[i],outphenos[i]) for i in range(len(outtypes)) ]

    to the end of the log
    containing tab separated filepaths and types which we now need to create in
    the history
    This code was written for the security branch and breaks if run in the 
    current trunk ! 
    """
    killme = string.punctuation + string.whitespace
    trantab = string.maketrans(killme,'_'*len(killme))
    job_name = param_dict.get( 'title', 'reArrayExpress Test' ).translate(trantab)
    aeid = param_dict.get('AEid','?')
    dbkey = param_dict.get('dbkey','hg18')
    base_dataset = out_data.items()[0][1]
    history = base_dataset.history
    if history == None:
        print "unknown history!"
        return
    hpath = out_data['htmlfile'].extra_files_path
    logfname = os.path.join(hpath,'%s.R.log' % job_name ) # determined in RRun!
    loglist = file(logfname,'r').readlines()
    logf = file(logfname,'a')
    newfiles = [x for x in loglist if x.split('\t')[0] == '### makenewgalaxy']
    newfiles = [x.strip().split('\t')[1:] for x in newfiles] # get rid of #makenewgalaxy
    logf.write('# read %s from %s' % (newfiles,logfname))
    # need pheno for each distinct name returned by arrayexpress
    for (file_path,file_type,pheno_path) in newfiles:
        # note, this gets created - pass extra args?
        #class DatasetInstance( object ):
        #    """A base class for all 'dataset instances', HDAs, LDAs, etc"""
        #    states = Dataset.states
        #    permitted_actions = Dataset.permitted_actions
        #    def __init__( self, id=None, hid=None, name=None, info=None, blurb=None, peek=None, extension=None, 
        #                  dbkey=None, metadata=None, history=None, dataset=None, deleted=False, designation=None,
        #                  parent_id=None, validation_errors=None, visible=True, create_dataset = False ):
        html = [galhtmlprefix % 'reArrayExpress',] # we need to attach the contents here
        newname = '%s.%s' % (aeid,file_type)
        html.append('<li><a href="%s">%s</a></li>' % (newname,newname))
        info = '%s, %s' % (job_name, aeid)
        pk = ''.join(file(pheno_path,'r').readlines())
        newdata = app.model.HistoryDatasetAssociation(extension=file_type,dbkey=dbkey,info=info,
            name=newname,peek = pk, create_dataset = True ) 
        # as noted in encode_import_code.py on which this was based :)
        # This import should become a library
        #newdata.metadata.base_name = aeid
        efp = newdata.extra_files_path
        try:
            os.makedirs(efp)
        except:
            pass
        phenoname = os.path.split(pheno_path)[-1] # name 
        newppath = os.path.join(efp,phenoname)
        shutil.copy(pheno_path,newppath) # save pheno for metadata
        html.append('<li><a href="%s">%s</a></li>' % (phenoname,phenoname))
        try:
            app.security_agent.set_dataset_permissions( newdata.dataset, base_dataset.dataset.groups )
        except:
            pass # old pre-security?
        app.model.context.flush()
        try:
            shutil.copyfile(file_path,os.path.join(efp,newname)) # copy into our extra files path 
            newdata.state = jobs.JOB_OK
        except:
            s = "The requested file %s is missing from the system." % file_path
            logf.write(s)
            logf.write('\n')
            logf.write('Trying to write to %s\n' % (newdata.file_name))
            newdata.info = s
            newdata.state = jobs.JOB_ERROR
        newdata.metadata.pheno_path = newppath  
        newdata.metadata.base_name = aeid      
        newdata.extension = file_type
        newdata.name = newname
        newdata.info = info
        newdata.set_peek()
        newdata.set_meta() # must set peek first
        logf.write('## saving %s as %s\n' % (newname, newdata.file_name))
        s = '# newdata %s peek = %s\n' % (newname,newdata.peek)
        logf.write(s)
        s = '# newdata %s metadata pheno_path = %s\n' % (newname,newdata.metadata.pheno_path)
        logf.write(s) 
        logf.write('\n')
        newdata.dbkey = dbkey
        newdata.set_size()
        html.append(galhtmlpostfix)
        f = open(newdata.get_file_name(),'w') # use inbuilt accessor
        f.write('\n'.join(html))
        f.write('\n')
        f.close()
        history.add_dataset( newdata )
        # newdata = app.datatypes_registry.change_datatype(newdata, file_type)
        app.model.context.flush()
    logf.write('post hook ends\n')
    logf.close()
    
        
            


