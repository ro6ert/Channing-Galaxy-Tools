# reMakeNew_code.py
# import this into an rexpression tool for a standard post hook
# to build any new rexpression datatypes found in the extra_files_path
# if no primary output html file content is found, creates new listing

from galaxy import datatypes, config, jobs
import logging,shutil,sys,string,os
import galaxy.util

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

gal_Log = logging.getLogger(__name__)

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
   res = [('no values found for phename %s' % phename,'?',False),]
   gal_Log.debug('##get_phevals called with phename=%s,rawphe=%s' % (phename,rawphe))
   for p in rawphe:
       if p[0] == phename:
           ares = []
           vlist = p[1] # (value,n) tuples
           for x in vlist:
                nocom = x[0].replace(',',oursep)
                ares.append(('%s: %s (%d)' % (phename,x[0],x[1]), nocom ,False))
                # encode commas in the value strings - these will be passed as params
                # so must be fixed on the way back in
                # note Sept09 that reGeoquery and reArrayExpress now
                # replace ;,: etc with - in all rows of the phenodata and remake the eset!
                # may need a fixer tool
           if len(ares) > 1:
           	res = ares
           break
   return res


def exec_after_process(app, inp_data, out_data, param_dict, tool, stdout, stderr):
    """    tricky. We need to do this for tools that return multiple new datasets or for
    generic tools that produce a bunch of things to be linked into a page
    If there are (eg) esets in the output, a new datatype should be created
    get output html extra_files_path file list as
    tab separated filepaths and types. We need links to those (or do we - > 2GB will kill the server
    when clicked as at Sept 2009, pointing to newly create datatypes in the history
    was in lots of code files - needed to be done once and right..
    Sacrifice ability to add application specific descriptive stuff to the names
    """
    firstPhepath = None # keep track of the first .pheno file we find
    re_exts = ['.eset','.affybatch','.malist']
    phe_ext = '.pheno'
    killme = string.punctuation + string.whitespace
    trantab = string.maketrans(killme,'_'*len(killme))
    job_name = param_dict.get( 'title', 'Job Title' ).translate(trantab)
    dbkey = param_dict.get('dbkey','hg18')
    base_dataset = out_data.items()[0][1]
    base_dataset.name = job_name
    history = base_dataset.history
    bdpath = base_dataset.file_name
    if history == None:
        gal_Log.error("## reMakeNew end job hook error on dataset %s - unknown history!" % (self.file_name))
        sys.exit(1)
    hpath = base_dataset.extra_files_path
    try:
         os.makedirs(hpath)
    except:
         pass
    try:
        newhtml = open(bdpath,'r').readlines()
        if len(newhtml) < 2:
            newhtml = [] # force recreate - too short - we want valid html
    except:
        newhtml = []
    flist = os.listdir(hpath)
    flist.sort()
    buildnewhtml = False
    if len(newhtml) == 0: # assume we have to write one
        buildnewhtml = True
        newhtml = [galhtmlprefix % job_name,] # we need to attach the contents here
        newhtml.append('Files created by this job:<br/>')
        newhtml.append('<ol>\n')
    gal_Log.info('## post hook buildnewhtml = %s' % buildnewhtml)
    makeme = []
    loglist = [] # use a new log if we don't find one
    logfname = os.path.join(hpath,'%s.log' % job_name)
    for f in flist:
        if buildnewhtml:
            newhtml.append('<li><a href="%s">%s</a></li>' % (f,f))
        fext = os.path.splitext(f)[-1]
        if fext in re_exts:
            makeme.append(f)
        if fext == '.log': # recycle the existing log
            loglist = file(os.path.join(hpath,f),'r').readlines() # recover log
            logfname = os.path.join(hpath,f)
    if buildnewhtml:
        newhtml.append('</ol></br>\n')
        newhtml.append(galhtmlpostfix)
        f = open(bdpath,'w') # use inbuilt accessor
        f.write('\n'.join(newhtml))
        f.write('\n')
        f.close()
    fpaths = [os.path.join(hpath,x) for x in makeme]
    flexts = [os.path.splitext(x)[-1][1:] for x in makeme] # use this metadata as datatype - without the leading .
    flnames = [os.path.splitext(x)[0] for x in makeme] # for pheno name
    fldict = dict(zip(flexts,range(len(flexts)))) # indices
    newfiles = [(fpaths[i],flexts[i],flnames[i]) for i in range(len(makeme))]
    gal_Log.info('## post hook newfiles = %s' % newfiles)
    logf = file(logfname,'w')
    if len(loglist) > 0:
        logf.write(''.join(loglist))
    logf.write('# found %s in the efp, and read %d rows from logfile %s' % (flist,len(loglist),logfname))
    # need pheno for each distinct name returned by arrayexpress
    for i, ftuple in enumerate(newfiles):
        (file_path,file_type,file_base) = ftuple
        alreadyExists = False
        if i > 0:
            newdata = app.model.HistoryDatasetAssociation(create_dataset = True, sa_session = app.model.context )
            hpath = newdata.extra_files_path
        else:
            alreadyExists = True # first one already exists - don't try to make
            newdata = base_dataset # is first one - already built
        # note, this gets created - pass extra args?
        #class DatasetInstance( object ):
        #    """A base class for all 'dataset instances', HDAs, LDAs, etc"""
        #    states = Dataset.states
        #    permitted_actions = Dataset.permitted_actions
        #    def __init__( self, id=None, hid=None, name=None, info=None, blurb=None, peek=None, extension=None,
        #                  dbkey=None, metadata=None, history=None, dataset=None, deleted=False, designation=None,
        #                  parent_id=None, validation_errors=None, visible=True, create_dataset = False ):
        newhtml = [galhtmlprefix % job_name,] # we need to attach the contents here
        newname = '%s.%s' % (file_base,file_type)
        newhtml.append('<li><a href="%s" type="application/octet-stream">%s</a></li>' % (newname,newname))
        info = '%s, %s' % (job_name, newname)
        phenoname = '%s%s' % (file_base,phe_ext)
        phenopath = os.path.join(hpath,phenoname) # should exist for all new rexpression datatypes?
        gal_Log.debug('### phenopath = %s' % phenopath)
        try:
            pk = ''.join(file(phenopath,'r').readlines())[:5]
            phefound = True
            if not firstPhepath:
                firstPhepath = phenopath # for 'ron' - later on
        except:
            if firstPhepath: # ok, let's recycle it
                pk = ''.join(file(firstPhepath,'r').readlines())[:5]
                phefound = True
                phenopath = firstPhepath # for the copy later
            else:
                phefound = False
                pk = 'No pheno file %s found for %s' % (phenopath,file_base)
                gal_Log.warning('## reMakeNew_code %s' % pk)
        # as noted in encode_import_code.py on which this was based :) This import should become a library
        newdata.extension=file_type # belt and braces
        newdata.dbkey=dbkey
        newdata.info=info
        newdata.name=newname
        newdata.peek=pk
        efp = newdata.extra_files_path
        newppath = os.path.join(efp,phenoname)
        newdata.metadata.pheno_path = newppath
        newdata.metadata.base_name = file_base
        if not alreadyExists:
            app.model.context.add( newdata )
            app.model.context.flush()
            app.security_agent.copy_dataset_permissions( base_dataset.dataset, newdata.dataset )
            history.add_dataset( newdata )
            app.model.context.add( history )
            app.model.context.flush()
            try:
                os.makedirs(efp)
            except:
                pass
            try:
                shutil.copy(phenopath,newppath) # save pheno for metadata
                newhtml.append('<li><a href="%s">%s</a></li>' % (phenoname,phenoname))
            except:
                gal_Log.error('## reMakeNew_code error - shutil.copy(%s,%s) failed' % (phenopath,newppath))
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
            logf.write('## saving %s as %s\n' % (newname, newdata.file_name))
            s = '# newdata %s peek = %s\n' % (newname,newdata.peek)
            logf.write(s)
            s = '# newdata %s metadata pheno_path = %s\n' % (newname,newdata.metadata.pheno_path)
            logf.write(s)
            logf.write('\n')
            newhtml.append(galhtmlpostfix)
            f = open(newdata.file_name,'w')
            f.write('\n'.join(newhtml))
            f.write('\n')
            f.close()
        newdata.set_meta() # must set peek first
        newdata = app.datatypes_registry.change_datatype(newdata, file_type)
    app.model.context.flush() # ? clean up
    logf.write('post hook ends\n')
    logf.close()
