#build list of available data
import os, sys, glob,shutil,time
import galaxy.util
from galaxy import datatypes, config, jobs 

mng = '### makenewgalaxy' # used by exec_after_process hook to parse out new file paths/names

repository = os.path.join(os.getcwd(),"tool-data/rg/library") # eeesh
#repository = "/home/rossl/rgalaxy/rgdata"


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

def timenow():
    """return current time as a string
    """
    return time.strftime('%d/%m/%Y %H:%M:%S', time.localtime(time.time()))

def get_available_file_types():
    import_files = []
    try:
        flist = glob.glob(os.path.join(repository,'*'))
        flist = [os.path.splitext(x)[0] for x in flist] # get unique filename
        flist = list(set(flist)) # remove dupes  
        flist.sort()
        for i, data in enumerate( flist ):
            import_files.append( (os.path.split(data)[-1], data, False) )
    except:
        pass
    if len(import_files) < 1:
        import_files.append(('No data available','None',True))
    return import_files

#return available datasets for build
def get_available_data( file_type_dir, build='hg18' ):
    #we need to allow switching of builds, and properly set in post run hook
    import_files = []
    try:
        flist = glob.glob(os.path.join(repository,file_type_dir,'*.ped'))  
        flist = [os.path.splitext(x)[0] for x in flist] # get unique filenames
        flist = list(set(flist)) # remove dupes  
        flist.sort()
        for i, fname in enumerate( flist ):
            import_files.append( ( os.path.split(fname)[-1], fname, False) ) # file, basename
    except:
        pass
    if len(import_files) < 1:
        import_files.append(('No data available for this build','None',True))
    return import_files


def exec_after_process(app, inp_data, out_data, param_dict, tool, stdout, stderr):
    """tricky. We wrote a list of basenames and the tdir to the end of the log that
    now need to be moved into pbed files in the history
    This code was written for the security branch 
    """
    dbkey = param_dict.get('dbkey','hg18')
    base_dataset = out_data.items()[0][1]
    history = base_dataset.history
    if history == None:
        print "unknown history!"
        return
    logpath = out_data['log'].file_name
    loglist = file(logpath,'r').readlines()
    newfiles = [x for x in loglist if x.split('\t')[0] == mng]
    newfiles = [x.strip().split('\t')[1:] for x in newfiles] # get rid of #makenewgalaxy
    d = os.getcwd()
    s = '### exec_after_process in %s. Newfiles = %s\n' % (d,newfiles)
    lf = file(logpath,'a')
    lf.write(s)
    lf.close()
    # parse out the encoded new datasets for galaxy
    for alist in newfiles:
        # note, this gets created - pass extra args?
        #class DatasetInstance( object ):
        #    """A base class for all 'dataset instances', HDAs, LDAs, etc"""
        #    states = Dataset.states
        #    permitted_actions = Dataset.permitted_actions
        #    def __init__( self, id=None, hid=None, name=None, info=None, blurb=None, peek=None, extension=None, 
        #                  dbkey=None, metadata=None, history=None, dataset=None, deleted=False, designation=None,
        #                  parent_id=None, validation_errors=None, visible=True, create_dataset = False ):
        fpath,bname = alist
        fpath = fpath.strip()
        newflist = glob.glob(os.path.join(fpath,'%s.*' % bname))
        newname = bname.strip()
        info = '%s converted by rgLpedPbed at %s' % (bname,timenow())
        newdata = app.model.HistoryDatasetAssociation(extension='pbed',dbkey=dbkey,info=info,
            name=newname, create_dataset = True )
        efp = newdata.extra_files_path
        try:
            os.makedirs(efp)
        except:
            pass
        
        # as noted in encode_import_code.py on which this was based :)
        # This import should become a library
        #newdata.metadata.base_name = geoid
        newdata.metadata.base_name = bname
        try:
            app.security_agent.set_dataset_permissions( newdata.dataset, base_dataset.dataset.groups )
        except:
            pass # old pre-security?
        app.model.context.flush()
        newdata.state = jobs.JOB_OK
        for f in newflist: # copy into efp or into the html file depending on extension
            n,e = os.path.splitext(f) # get extension
            if e.lower() == '.html': # contents of the base html file
                try:
                    shutil.copyfile(f,newdata.file_name) 
                    #os.unlink(f) # get rid of source
                except:
                    s = "The requested file %s is missing from %s." % (f,file_path)
                    lf = file(logpath,'a')
                    lf.write(s)
                    lf.write('\n')
                    lf.write('Trying to write %s to %s\n' % (f,newdata.file_name))
                    lf.close()
                    newdata.info = s
                    newdata.state = jobs.JOB_ERROR
            else: # one of the basename files in the tempdir
                 newf = os.path.split(f)[-1] # get name
                 newfpath = os.path.join(efp,newf) # new name for move
                 try:
                    shutil.move(f,newfpath) # move into extra files path 
                 except:
                    s = "The requested file %s is missing from %s." % (f,file_path)
                    lf = file(logpath,'a')
                    lf.write(s)
                    lf.write('\n')
                    lf.write('Trying to write %s to %s\n' % (f,efp))
                    lf.close()
                    newdata.info = s
                    newdata.state = jobs.JOB_ERROR
               
        newdata.dbkey = dbkey
        newdata.set_peek()
        newdata.set_meta() # must set peek first
        newdata.set_size()
        history.add_dataset( newdata )
        app.model.context.flush()
        doned = {}
##        for alist in newfiles: # clean up
##            fpath,bname = alist
##            if not doned.get(fpath,None):
##                try:
##                    shutil.rmtree(fpath)
##                    doned[fpath] = fpath
##                except:
##                    pass
##            


