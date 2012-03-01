#TODO: Set dbkey to proper UCSC build, if known
# Greg Von Kuster
import urllib

from galaxy import datatypes, config
import tempfile, shutil

def exec_before_job( app, inp_data, out_data, param_dict, tool=None):
    """Sets the name of the data"""
    data_name = param_dict.get( 'name', 'Biomart query' )
    data_type = param_dict.get( 'type', 'text' )
    if data_type == 'text': data_type='interval' #All data from biomart is TSV, assume interval
    name, data = out_data.items()[0]
    data = app.datatypes_registry.change_datatype(data, data_type)
    data.name = data_name
    out_data[name] = data

def exec_after_process(app, inp_data, out_data, param_dict, tool=None, stdout=None, stderr=None):
    """Verifies the data after the run"""

    URL = param_dict.get( 'URL', None )
    URL = URL + '&_export=1&GALAXY_URL=0'
    if not URL:
        raise Exception('Datasource has not sent back a URL parameter')

    CHUNK_SIZE = 2**20 # 1Mb 
    MAX_SIZE   = CHUNK_SIZE * 100
    
    try:
        # damn you stupid sanitizer! - URL is now in the NEVER_SANITIZE list in util.py
        #URL = URL.replace('martX', 'mart&')
        #URL = URL.replace('0X_', '0&_')
        page = urllib.urlopen(URL)
    except Exception, exc:
        raise Exception('Problems connecting to %s (%s)' % (URL, exc) )

    name, data = out_data.items()[0]
    
    fp = open(data.file_name, 'wb')
    size = 0
    while 1:
        chunk = page.read(CHUNK_SIZE)
        if not chunk:
            break
        if size > MAX_SIZE:
            raise Exception('----- maximum datasize exceeded ---')
        size += len(chunk)
        fp.write(chunk)

    fp.close()
    #Set meta data, format file to be valid interval type
    if isinstance(data.datatype, datatypes.interval.Interval):
        data.set_meta(first_line_is_header=True)
        #check for missing meta data, if all there, comment first line and process file
        if not data.missing_meta():
            line_ctr = -1
            temp = tempfile.NamedTemporaryFile('w')
            temp_filename = temp.name
            temp.close()
            temp = open(temp_filename,'w')
            chromCol = int(data.metadata.chromCol) - 1
            startCol = int(data.metadata.startCol) - 1
            strandCol = int(data.metadata.strandCol) - 1
            
            
            for line in open(data.file_name, 'r'):
                line_ctr += 1
                
                #First line is a non-commented header line, lets comment it out here
                if line_ctr == 0:
                    temp.write("#%s" % line)
                    continue
                
                fields = line.strip().split('\t')
                
                #If chrom col is an int, make it chrInt
                try:
                    int(fields[chromCol])
                    fields[chromCol] = "chr%s" % fields[chromCol]
                except:
                    try:
                        if fields[chromCol].upper()== "X" or fields[chromCol].upper()== "Y":
                            fields[chromCol] = "chr%s" % fields[chromCol].upper()
                    except:
                        pass
                    
                #change to BED coordinate system
                try:
                    fields[startCol] = str(int(fields[startCol]) - 1)
                except:
                    pass
                
                #set strand to +/-, instead of +1/-1
                try:
                    if strandCol > 0:
                        if int(fields[strandCol]) > 0:
                            fields[strandCol] = "+"
                        else:
                            fields[strandCol] = "-"
                except:
                    pass
                
                temp.write("%s\n" % '\t'.join(fields))
            
            temp.close()
            shutil.move(temp_filename,data.file_name)
            
        else:
            data = app.datatypes_registry.change_datatype(data, 'tabular')
    data.set_peek()
    app.model.context.flush()
