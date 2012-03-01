# before running the qc, need to rename various output files
import os,string
from galaxy import datatypes

def get_phecols(phef):
   """return column names """
   head = file(phef.file_name,'r').next()
   c = head.strip().split()
   res = [(cname,cname,False) for cname in c]
   if len(res) > 0:
      x,y,z = res[0] # 0,1 = fid,iid
      res[0] = (x,y,True) # set second selected
   else:
      res = [('','no columns found',False),]
   return res

def exec_after_process(app, inp_data, out_data, param_dict, tool, stdout, stderr):
    """Sets the name of the data"""
    trantab = string.maketrans(string.punctuation,'_'*len(string.punctuation))
    job_name = param_dict.get( 'title1', 'MultiSort' ).translate(trantab)
    newname = '%s_MultiSort.xls' % job_name
    n = 'out_file1'
    data = out_data['out_file1']
    data = app.datatypes_registry.change_datatype(data, 'tabular')
    data.name = newname
    out_data[n] = data

