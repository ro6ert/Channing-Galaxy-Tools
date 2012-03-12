#build list of available data
import string, os, sys, glob, shutil
import galaxy.util
from galaxy import datatypes, config, jobs 


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

def exec_after_process(app, inp_data, out_data, param_dict, tool, stdout, stderr):
    """tricky. We wrote     
    s = paste(ngs,'\t',outf_maList,'\tmalist\t',outf_pheno)
    cat(s)
    to the end of the log
    containing tab separated filepaths and types which we now need to create in
    the history
    This code was written for the security branch 
    """
    mng = '### makenewgalaxy'
    killme = string.punctuation + string.whitespace
    trantab = string.maketrans(killme,'_'*len(killme))
    job_name = param_dict.get( 'title', 'reGeoQuery Test' ).translate(trantab)
    geoid = param_dict.get('gid','?')
    dbkey = param_dict.get('dbkey','hg18')
    base_dataset = out_data.items()[0][1]
    history = base_dataset.history
    if history == None:
        print "unknown history!"
        return
    logpath = out_data['logmeta'].file_name
    loglist = file(logpath,'r').readlines()
    newfiles = [x for x in loglist if x.split('\t')[0] == mng]
    # parse out the encoded new datasets for galaxy
    newfiles = [x.strip().split('\t')[1:] for x in newfiles] # get rid of #makenewgalaxy
    phenonamesmade = {}
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
        file_path = file_path.strip() # eg /tmp/GSE6536dJbdcI/GSE6536_1.eset 
        file_type = file_type.strip()
        newname = os.path.split(file_path)[-1]
        #newname = '%s.%s' % (geoid,file_type)
        info = '%s, %s, %s' % (job_name, geoid, newname)
        pheno_path = pheno_path.strip()
        if pheno_path <> 'None': # an rexpression datatype
            pk = file(pheno_path,'r').read()
            newdata = app.model.HistoryDatasetAssociation(extension=file_type,dbkey=dbkey,info=info,
                name=newname,peek = pk[:500], create_dataset = True ) 
            # as noted in encode_import_code.py on which this was based :)
            # This import should become a library
            #newdata.metadata.base_name = geoid
            efp = newdata.extra_files_path
            try:
                os.makedirs(efp)
            except:
                pass
            phenoname = os.path.split(pheno_path)[-1] # name 
            newppath = os.path.join(efp,phenoname)
            shutil.copy(pheno_path,newppath) # save pheno for metadata
            newdata.metadata.pheno_path = newppath  
            newdata.metadata.base_name = geoid
        else: # is (eg) a tabular pheno file
            newdata = app.model.HistoryDatasetAssociation(extension=file_type,dbkey=dbkey,info=info,
                name=newname, create_dataset = True )         
        newdata.extension = file_type
        newdata.name = newname
        newdata.info = info
        try:
            app.security_agent.set_dataset_permissions( newdata.dataset, base_dataset.dataset.groups )
        except:
            pass # old pre-security?
        app.model.context.flush()
        try:
            shutil.copyfile(file_path,newdata.file_name) 
            newdata.state = jobs.JOB_OK
        except:
            s = "The requested file %s is missing from the system." % file_path
            lf = file(logpath,'a')
            lf.write(s)
            lf.write('\n')
            lf.write('Trying to write to %s\n' % (newdata.file_name))
            lf.close()
            newdata.info = s
            newdata.state = jobs.JOB_ERROR
        newdata.dbkey = dbkey
        newdata.set_peek()
        newdata.set_meta() # must set peek first
        lf = file(logpath,'a')
        lf.write('## saving %s as %s\n' % (newname, newdata.file_name))
        s = '# newdata %s peek = %s\n' % (newname,newdata.peek)
        lf.write(s)
        s = '# newdata %s metadata pheno_path = %s\n' % (newname,newdata.metadata.pheno_path)
        lf.write(s) 
        lf.write('\n')
        lf.close()
        newdata.set_size()
        history.add_dataset( newdata )
        app.model.context.flush()
