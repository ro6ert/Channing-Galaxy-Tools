#!/usr/bin/env python2.4
"""
june 2009 - changes to deal with cnv variants in illumina data
can't use array conversion without removing multi-allelic markers
"""

import optparse
import os
import gzip
import commands
from sets import Set
from array import array

try:
    import MySQLdb
    USE_DB_CACHE = True
except ImportError:
    print 'WARNING: Cannot import one or more modules required for database caching.'
    USE_DB_CACHE = False
    
VERBOSE = True

MISSING_ALLELES = Set(['N', '0', '.', '-', 0])

BYTES_PER_GENOTYPE = 4

AUTOSOMES = set(range(1, 23) + [str(c) for c in range(1, 23)])

class DuplicateMarkerInMapFile(Exception): pass
class MapLineHasThreeFields(Exception): pass

MAP_LINE_EXCEPTION_TEXT = """
One or more lines in the *.map file has only three fields.
The line was:

%s

If you are running rgGRR through EPMP, this is usually a
sign that you are using an old version of the map file.
You can correct the problem by re-running Subject QC.  If
you have already tried this, please contact the developers,
or file a bug.
"""

class Ped:
    def __init__(self, path=None, doParse=False):
        self.path = path
        self.name = None
        self.markers = Set()
        self._ordered_markers = []
        self._ordered_subjects = []
        self._subject_keys_by_line = {}
        self.subjects = Set()
        self.samples = Set()
        self._subject_info = {}
        self.nGenotypes = 0
        self._genotypes = {}
        #self._INDEXES_BY_SUBJ = {}
        #self._INDEXES_BY_LINE = []
        self._fh = None
        self.badmarkers = [] # multiallelic or indel marker index list
        if os.path.exists(path):
            self.preparse()
            if doParse:
                self.bad = self.parse()
        else:
            raise 'NO file %s exists' % (path), str(os.listdir('.'))

    def preparse(self):
        """ Pre-parse the marker and subject information for this ped-file
        """

        ### md5sum 
        ### md5sum = commands.getoutput('md5sum %s' % (self.path)).split()[0].strip()
        md5sum = 'md5sum here'
        
        ### First the map file
        if VERBOSE:
            print 'pedfile.PED: "%s" (%s)' % (os.path.basename(self.path), md5sum)

        ### process the companion map file
        map_path = self.path.replace('.ped', '.map')
        if os.path.exists(map_path):
            self.mapf = Map(map_path)
        else:
            raise 'NoMap'
        nMarkers = len(self.mapf)
        if VERBOSE:
            print '\tlen(rslist)=%s' % (nMarkers)

        bytesBeforeLine = 0
        bytesForRestOfLine = nMarkers * BYTES_PER_GENOTYPE
        pedf = self._fh = open(self.path, 'r')
        for lnum, line in enumerate(pedf):
            #fid, oiid, odid, omid, sex, phe, genos = line.strip().split(None, 6)
            Lline = line.strip().split()
            fid, oiid, odid, omid, sex, phe = Lline[:6]
            genos = Lline[6:]
            header = (fid, oiid, odid, omid, sex, phe)
            iid = oiid.split('.')[0]
            did = odid.split('.')[0]
            mid = omid.split('.')[0]
            
            key = (fid, iid)
            okey = (fid, oiid)
            
            self.subjects.add(okey)
            self._ordered_subjects.append(okey)
            self._subject_info[okey] = (fid, iid, did, mid, sex, phe, oiid)
            self._subject_keys_by_line[lnum] = okey
            ### print '\tSUBJECT: %s' % (str(header))
        print 'MARKERS:   %s' % (len(self.mapf))
        print 'SUBJECTS:  %s' % (len(self.subjects))

    def getGenotypeByLineAndIndex(self, line, index):
        """ given a row, retrieve the index'th genotype
        """
        subject_key = self._subject_keys_by_line[line]
        gindex = index * 2
        genotypes = self._genotypes.get(subject_key,None)
        genotype = ['0','0']
        if (genotypes <> None):
            if (len(genotypes) > (gindex+2)):
                genotype = genotypes[gindex:(gindex+2)]
            else:
                print '### len=%d for genotypes %s wanting %d' % (len(genotypes),','.join(subject_key),gindex+2)
        #else:
        #    print '### No genos for for subject %s' % (','.join(subject_key))
        return genotype

    def getSubjectInfo(self, fid, oiid):
        """
        """
        return self._subject_info[(fid, oiid)]

    def getSubjectInfoByLine(self, line):
        """
        """
        return self._subject_info[self._ordered_subjects[line]]
        
    def parse(self, path=None):
        """ Parse a given file -- this needs to be memory-efficient so that large
            files can be parsed (~1 million markers on ~5000 subjects?)
            need to deal with indels/cnvs - initially by setting to null
            probably possible to maintain a dict of alleles and fake as (eg g/t) snp
            but for now - set all these to missing
        """
        print '### parse'
        if path:
            self.path = path

        if self.path.endswith('.gz'):
            pfile = gzip.open(self.path, 'r')
        else:
            pfile = open(self.path, 'r')

        mapf = self.mapf
        badcols = {}
        nLines = len(self.subjects)
        for lnum, line in enumerate(pfile):
            line = line.strip()
            if not line:
                continue
            Lline = line.split()
            fid, oiid, odid, omid, sex, phe = Lline[:6]
            genos = Lline[6:]
            for i,g in enumerate(genos): # hack to set indel/cnv/multis to missing
                if badcols.get(i/2,None) or (g == '-') or (len(g) > 1): # indel or multiallelic
                    bc = i/2 # snp# - integer div - make sure we get both of each pair?
                    genos[2*bc] = '0' # set to missing for now
                    genos[2*bc + 1] = '0' # both alleles for this snp
                    if not badcols.get(bc,None): # first time
                        badcols.setdefault(bc,1) # make sure all are set to missing
                        print '## bad col %d = %s' % (bc,self.mapf.get_index_marker(bc))
            k = self._subject_keys_by_line[lnum]
            self.nGenotypes += len(genos)/2
            try:
                self._genotypes[k] = array('c', genos)
            except:
                print '### bad array cast for %s' % k
                print '### genos = %s' % genos[:500]
        bk = badcols.keys()
        if len(bk) > 0:
            bk.sort()
        self.badmarkers = bk # save for processing
        print '#GENOTYPES: %d' % (self.nGenotypes)
        print '#Bad genotypes: %d' % (len(bk))

    def badmarker_rep(self):
        res = []
        if len(self.badmarkers) > 0:
            res = ['%s\t%d' % (self.mapf.get_index_marker(x),x) for x in bk]
        return res
        

    def autosomal_indices(self):
        """ Return the indices of markers in this ped/map that are autosomal.
            This is used by rgGRR so that it can select a random set of markers
            from the autosomes (sex chroms screw up the plot)
        """
        return self.mapf._autosomal_indices

class Map:
    def __init__(self, path=None):
        self.path = path
        self.markers = {}
        self._ordered_markers = []
        self._autosomal_indices = set()
        if os.path.exists(path):
            self.parse()

    def __len__(self):
        return len(self.markers)

    def get_index_marker(self,index):
        """ return marker name for this index
        """
        return self._ordered_markers[index]

    def parse(self):
        ndupes = 0
        print 'MAP: %s' % (self.path)
        warned = False
        if self.path.endswith('.gz'):
            fh = gzip.open(self.path, 'r')
        else:
            fh = open(self.path, 'r')
        for i, line in enumerate(fh):
            line = line.strip()
            if not line:
                continue

            fields = line.split()
            if len(fields) == 1:
                if not warned:
                    print '\tWARNING: This looks like one of those broken map files with only snp names'
                    warned = True
                name = fields[0]
                chrom = -1
                genpos = 0
                abspos = 0
            elif len(fields) == 3:
                raise MapLineHasThreeFields(MAP_LINE_EXCEPTION_TEXT % (str(line)))
            else:
                chrom, name, genpos, abspos = fields
                if chrom in AUTOSOMES:
                    self._autosomal_indices.add(i)
            if name in self.markers:
                #raise DuplicateMarkerInMapFile('Marker %s was found twice in map file %s' % (name, self.path))
                ndupes += 1
                name = '%s_%d' % (name,ndupes) # make a unique name
            self.markers[name] = (chrom, name, genpos, abspos)
            self._ordered_markers.append(name)
        fh.close()
        
            
if __name__ == '__main__':
    """ Run as a command-line, this script should get one or more arguments,
        each one a ped file to be parsed with the PedParser (unit tests?)
    """
    op = optparse.OptionParser()
    opts, args = op.parse_args()
    for arg in args:
        ped = Ped(arg,doParse=True)


