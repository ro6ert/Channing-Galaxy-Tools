from galaxy import app
import galaxy.util

librepos = '/usr/local/galaxy/data/rg'
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

def getcols(fname="/usr/local/galaxy/data/camp2007/camp2007.xls"):
   """return column names other than chr offset as a select list"""
   head = open(fname,'r').next()
   c = head.strip().split()
   res = [(cname,'%d' % n,False) for n,cname in enumerate(c)]
   for i in range(4):
     x,y,z = res[i]
     res[i] = (x,y,True) # set first couple as selected
   return res


def exec_after_process(app, inp_data, out_data, param_dict, tool, stdout, stderr):
    """Sets the name of the data"""
    data_name = param_dict.get( 'title', 'My Camp2007 results' )
    name, data = out_data.items()[0]
    data = app.datatypes_registry.change_datatype(data, 'fbat')
    data.name = data_name
    out_data[name] = data

