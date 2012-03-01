import os, sys, glob, time

lookup = {}
lookup['lped'] = 'Linkage (separate ped/map) format'
lookup['eigenstrat'] = 'Eigenstrat format'
lookup['fped'] = 'FBat format (Marker names in header row)'
lookup['pbed'] = 'Plink compressed binary format'
lookup['snpmatrix'] = 'snpMatrix R .rda binary format'
lookup['pphe'] = 'Plink phenotype file'
lookup['fphe'] = 'Fbat phenotype file'
lookup['eset'] = 'Rexpression eSet file'
lookup['affybatch'] = 'Rexpression AffyBatch file'

#GALAXY_DATA_INDEX_DIR = app.config.tool_data_path    

def timenow():
    """return current time as a string
    """
    return time.strftime('%d/%m/%Y %H:%M:%S', time.localtime(time.time()))

def get_available_file_types(REPOSITORY=None):
    import_files = []
    try:
        flist = glob.glob(os.path.join(REPOSITORY,'*'))
        flist = [os.path.splitext(x)[0] for x in flist] # get unique dirnames
        flist = list(set(flist)) # remove dupes
        flist.sort()
        for i, data in enumerate( flist ):
            d = os.path.split(data)[-1]
            s = lookup.get(d,d) # translate or as is
            if s <> 'lped':
                import_files.append( (s, data, False) )
            else:
                import_files.append( (s, data, True) )
    except:
	pass
    if len(import_files) < 1:
        import_files.append(('No system library files available','None',False))
    return import_files

#return available datasets for build
def get_available_data( REPOSITORY, file_ext, build='hg18' ):
    #we need to allow switching of builds, and properly set in post run hook
    import_files = []
    try:
        dname = os.path.join(REPOSITORY,file_ext,'*.*')
        flist = glob.glob(dname)
        flist = [os.path.splitext(x)[0] for x in flist] # get unique filenames
        flist = list(set(flist)) # remove dupes
        flist.sort()
        for i, data in enumerate( flist ):
            import_files.append( (os.path.split(data)[-1], os.path.split(data)[-1], False) )
    except:
	pass
    if len(import_files) < 1:
        import_files.append(('No data available for this build','None',True))
    return import_files



# Create link to files here
def exec_after_process(app, inp_data, out_data, param_dict, tool, stdout, stderr):
    name,data = out_data.items()[0]
    ftd = param_dict['file_type_dir']
    newdt = os.path.split(ftd)[-1]
    data.change_datatype(newdt)
    data.file_name = data.file_name
    base_name = param_dict['base_name']
    info = '%s created by rgenetics_import_code.py at %s' % (base_name,timenow())
    data.metadata.base_name = base_name
    data.name = '%s.%s' % (base_name,newdt)
    data.info = info
    app.model.context.flush()
 


