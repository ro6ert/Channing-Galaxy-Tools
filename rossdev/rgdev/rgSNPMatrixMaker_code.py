from galaxy import app
import os, string

def exec_after_process(app, inp_data, out_data, param_dict, tool, stdout, stderr):
    for name,data in out_data.items():
        if name == 'outfile1':
            ftype = 'snpmatrix'
            base_name = param_dict.get( 'title', 'My snpMatrix' ).encode()
            ss = '%s%s' % (string.punctuation,string.whitespace)
            ptran =  string.maketrans(ss,'_'*len(ss))
            base_name = base_name.translate(ptran)
            data.change_datatype(ftype)
            data.file_name = data.file_name
            data.metadata.base_name = base_name
            data.name = base_name
        elif name == 'logpath':
            data.name = '%s.log' % base_name
    app.model.context.flush()

