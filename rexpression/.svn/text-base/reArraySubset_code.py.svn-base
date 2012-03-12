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

oursep = '=+comma+=' # shouldn't be in any data?



def get_phecols(eset):
   """return phenotype column names for an rexpression phenodata as a dynamic
   select list for a Galaxy form - user selects one and then gets to select which values
   from the right get_phevals select.
   There is some convoluted logic in genetics.py to
   generate the pheCols metadata - do it once rather than every time we need it
   A column with only 1 value doesn't change, so is not interesting for
   analysis. A column with a different value in every row is equivalent to a unique
   identifier so is also not interesting for anova or limma analysis so both these
   are removed after the concordance (count of unique terms) is constructed for each
   column. Then a complication - each remaining pair of columns is tested for
   redundancy - if two columns are always paired, then only one is needed"""
   rawphe = eset.metadata.pheCols
   t = type(rawphe)
   if t == type(u'xy') or t == type('xy') : # make a list
        rawphe = list(eval(rawphe))   
   # now res is [('foo',[('a',3),('b',11),..]),..]
   # eg [[u'tissue', [[u'spleen', 14], [u'thymus', 14]]],
   # WAS metadata.pheCols becomes [('bar;22, zot;113','foo'), ...]
   # WAS res = [('Col:%s has %s' % (x[1],x[0]),x[1],False) for x in rawphe]
   res = []   
   for p in rawphe:
       if len(p) == 2:
         (cname,cvals) = p
         vs = ['%s (%d)' % (x[0],x[1]) for x in cvals]
         s1 = '%s contains %s' % (cname,'; '.join(vs))
         res.append((s1,cname,False))
       else:
         gal_Log.warning('### wierd - p = %s' % p) 
   if len(res) >= 1:
      x,y,z = res[0] # 0,1 = fid,iid
      res[0] = (x,y,True) # set one selected
   else:
      res = [('','no phenotype columns found',False),]
   return res

def get_phevals(eset,phename):
   """return phenotype values and counts for an rexpression phenotype column
   as a dynamic select list"""
   # now res is [['foo',[['a',3],['b',11],..]],..]
   # eg [[u'tissue', [[u'spleen', 14], [u'thymus', 14]]],
   rawphe = eset.metadata.pheCols
   t = type(rawphe)
   if t == type(u'xy') or t == type('xy') : # make a list
        rawphe = list(eval(rawphe))   
   res = [('?','no values found for phename %s' % phename,False),]
   for p in rawphe:
       if p[0] == phename:
           ares = []
           vlist = p[1] # (value,n) tuples
           for x in vlist:
               nocom = x[0].replace(',',oursep)
               ares.append('%s: %s (%d)' % (phename,x[0],x[1]), nocom ,False)
           # encode commas in the value strings - these will be passed as params so must be fixed on the
           # way back in 
           break
       if len(ares) > 1:
            res = ares
   return res



def exec_after_process(app, inp_data, out_data, param_dict, tool, stdout, stderr):
    """tricky. We wrote     
    outlist = ['%s\t%s\t%s' % (mng,newfiles[i],newnames[i]) for i in range(length(newfiles))]
    to the end of the log
    containing tab separated filepaths and types which we now need to create in
    the history
    This code was written for the security branch 
    """
    mng = '### makenewgalaxy'
    killme = string.punctuation + string.whitespace
    trantab = string.maketrans(killme,'_'*len(killme))
    job_name = param_dict.get( 'title', 'reGeoQuery Test' ).translate(trantab)
    iphenofile = inp_data.items()[0][1].metadata.pheno_path
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
    for (file_path,file_name,file_type) in newfiles:
        # note, this gets created - pass extra args?
        #class DatasetInstance( object ):
        #    """A base class for all 'dataset instances', HDAs, LDAs, etc"""
        #    states = Dataset.states
        #    permitted_actions = Dataset.permitted_actions
        #    def __init__( self, id=None, hid=None, name=None, info=None, blurb=None, peek=None, extension=None, 
        #                  dbkey=None, metadata=None, history=None, dataset=None, deleted=False, designation=None,
        #                  parent_id=None, validation_errors=None, visible=True, create_dataset = False ):
        file_path = file_path.strip()
        newname = file_name.strip()
        info = '%s, %s' % (job_name, newname)
        pk = file(iphenofile,'r').read()
        newdata = app.model.HistoryDatasetAssociation(extension=file_type,dbkey=dbkey,info=info,
            name=newname,peek = pk, create_dataset = True ) 
        # as noted in encode_import_code.py on which this was based :)
        # This import should become a library
        #newdata.metadata.base_name = geoid
        efp = newdata.extra_files_path
        try:
            os.makedirs(efp)
        except:
            pass
        phenoname = os.path.split(iphenofile)[-1] # name 
        iid = os.path.splitext(phenoname)[0]
        newppath = os.path.join(efp,phenoname)
        shutil.copy(iphenofile,newppath) # save pheno for metadata
        newdata.metadata.pheno_path = newppath  
        newdata.metadata.base_name = iid
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
