from galaxy import app
import galaxy.util,string,os,glob

librepos = '/usr/local/galaxy/data/rg/library'
myrepos = '/home/rerla/galaxy'
marchinirepos = '/usr/local/galaxy/data/rg/snptest'

#Provides Upload tool with access to list of available builds

builds = []
#Read build names and keys from galaxy.util
for dbkey, build_name in galaxy.util.dbnames:
    builds.append((build_name,dbkey,False))

#Return available builds
def get_available_builds(defval='hg18'):
    for i,x in enumerate(builds):
        if x[1] == defval:
           x = list(x)
           x[2] = True
           builds[i] = tuple(x)
    return builds

#return available datasets for build
def get_available_data( file_type_dir, build='hg18' ):
    #we need to allow switching of builds, and properly set in post run hook
    import_files = [('Use the History file chosen below','',True),]
    flist = glob.glob(os.path.join(librepos,file_type_dir,'*.ped'))  
    flist = [os.path.splitext(x)[0] for x in flist] # get unique filenames
    flist = list(set(flist)) # remove dupes  
    flist.sort()
    for i, data in enumerate( flist ):
        #import_files.append( (os.path.split(data)[-1], os.path.split(data)[-1], False) )
        import_files.append((data,data, False) )
    if len(import_files) < 1:
        import_files.append(('No %s data available - please choose a History file below'
                             % file_type_dir,'',True))
    return import_files

def get_available_outexts(uid):
    userId = uid
    print 'userId=',userId
    flib = os.path.join(librepos,userId,'fped')
    plib = os.path.join(librepos,userId,'lped')
    e = [('Fbat style (header row has marker names)',flib,False),('Linkage style (separate .map file)',plib,True)]
    return e

def getcols(fname="/usr/local/galaxy/data/camp2007/camp2007.xls"):
   """return column names other than chr offset as a select list"""
   head = open(fname,'r').next()
   c = head.strip().split()
   res = [(cname,'%d' % n,False) for n,cname in enumerate(c)]
   for i in range(4):
     x,y,z = res[i]
     res[i] = (x,y,True) # set first couple as selected
   return res

# Create link to files here
def exec_after_process(app, inp_data, out_data, param_dict, tool, stdout, stderr):
    base_name = param_dict.get( 'title', 'My Camp2007 results' )
    ptran =  string.maketrans(string.punctuation,'_'*len(string.punctuation))
    base_name = base_name.translate(ptran)
    name,data = out_data.items()[0]
    # uid = param_dict['userId']
    uid = '%d' % (data.history.user.id)
    # /usr/local/galaxy/data/rg/1/lped/
    ftd = os.path.join('/usr/local/galaxy/data/rg',uid,'lped')
    data.change_datatype('lped')
    data.file_name = data.file_name
    data.metadata.base_name = base_name
    data.name = base_name
    data.flush()
    data.dataset_file.extra_files_path = ftd
    data.dataset_file.readonly = True
    data.flush()
    app.model.flush()

