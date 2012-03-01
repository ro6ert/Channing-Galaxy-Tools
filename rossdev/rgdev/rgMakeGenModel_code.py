# before running the qc, need to rename various output files
import os,string
from galaxy import datatypes



def exec_after_process(app, inp_data, out_data, param_dict, tool, stdout, stderr):
    """Sets the name of the data"""
    trantab = string.maketrans(string.punctuation,'_'*len(string.punctuation))
    job_name = param_dict.get( 'title', 'GeneticModel' ).encode().translate(trantab)
    outxls = ['xls','%s.xls' % job_name]
    logtxt = ['txt','%s.log' % job_name]
    lookup={}
    lookup['output1'] = outxls
    lookup['logfile'] = logtxt
    for name in lookup.keys():
        data = out_data[name]
        data_type,newname = lookup[name]
        data = app.datatypes_registry.change_datatype(data, data_type)
        data.name = newname
        out_data[name] = data

