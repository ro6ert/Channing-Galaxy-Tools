"""
"""
from galaxy import datatypes,model
import sys,time,string

def timenow():
    """return current time as a string
    """
    return time.strftime('%d/%m/%Y %H:%M:%S', time.localtime(time.time()))


def exec_after_process(app, inp_data, out_data, param_dict, tool, stdout, stderr):
    """Change data file names

    """
    killme = string.punctuation + string.whitespace
    trantab = string.maketrans(killme,'_'*len(killme))
    job_name = param_dict.get( 'out_prefix', 'PicardHsMetrics' )
    title = job_name.translate(trantab)
    html = ['html','%s.html' % title]
    lookup={}
    lookup['html_file'] = html
    info = '%s Picard CalculateHsMetrics report created by rgPicardHsMetrics at %s' % (title,timenow())
    for aname in lookup.keys():
       data = out_data[aname]
       data_type,newname = lookup[aname]
       data = app.datatypes_registry.change_datatype(data, data_type)
       data.name = newname
       data.info = info
       out_data[aname] = data
    app.model.context.flush()



