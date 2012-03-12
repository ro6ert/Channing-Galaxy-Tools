#build list of available data
import string, os, sys, glob, shutil,logging
import galaxy.util
from galaxy import datatypes, config, jobs 


repository = "/usr/local/galaxy/data/rg/library"
gal_Log = logging.getLogger(__name__)

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

oursep = '=+comma+=' # shouldn't be in any data


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
            for i,x in enumerate(vlist):
                nocom = x[0].replace(',',oursep)
                ares.append(('%s: %s (%d)' % (phename,x[0],x[1]), nocom ,(i<2))) # set first 2 selected as default
            # encode commas in the value strings - these will be passed as params
            # so must be fixed on the way back in
            # note Sept09 that reGeoquery and reArrayExpress now
            # replace ;,: etc with - in all rows of the phenodata and remake the eset!
            # may need a fixer tool
            if len(ares) >= 1:
                res = ares
            break

    return res



def exec_after_process(app, inp_data, out_data, param_dict, tool, stdout, stderr):
    """Sets the name of the html file to the title"""
    trantab = string.maketrans(string.punctuation,'_'*len(string.punctuation))
    job_name = param_dict.get( 'title', 'ArrayTools' ).translate(trantab)
    data = out_data['outhtml']
    newname = '%s.html' % job_name
    data.name = newname
    out_data['outhtml'] = data
    app.model.context.flush()

