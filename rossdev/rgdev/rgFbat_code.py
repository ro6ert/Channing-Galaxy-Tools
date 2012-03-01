# before running the qc, need to rename various output files


from galaxy import datatypes,model
import os, string

def get_phecols(phef):
   """return column names """
   res = []
   if phef <> None:
       phepath = phef.extra_files_path
       phename = phef.name
       #phe = os.path.join(phepath,'%s.fphe' % (phename))
       head = open(phename,'r').next()
       c = head.strip().split()
       res = [(cname,cname,False) for cname in c]
   if len(res) > 0:
      x,y,z = res[0] # 0,1 = fid,iid
      res[0] = (x,y,True) # set second selected
   else:
      res = [('','no phenotype columns found',False),]
   return res


"""
      <command interpreter="python2.4">
        fbat.py $pedf $phef "$title" "$phecols" $out_file1 $logfile $smendelfile $mmendelfile
    </command> 
"""

def exec_after_process(app, inp_data, out_data, param_dict, tool, stdout, stderr):
    """Sets the name of the data"""
    trantab = string.maketrans(string.punctuation,'_'*len(string.punctuation))
    job_name = param_dict.get( 'title', 'Fbat' ).translate(trantab)
    dname = 'out_file1'
    data = out_data[dname]
    data_type = 'gg'
    newname = '%s_fbat.gg' % job_name
    data = app.datatypes_registry.change_datatype(data, data_type)
    data.name = newname
    out_data[dname] = data

