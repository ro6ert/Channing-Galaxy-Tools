""" Utility methods to be used in various lab projects
"""
import os,time
from sets import Set
import commands

### A list of IUPAC ambiguity codes
IUPAC = {
     'AC':'M',
     'AG':'R',
     'AT':'W',
     'CG':'S',
     'CT':'Y',
     'GT':'K',
     'ACG':'V',
     'ACT':'H',
     'AGT':'D',
     'CGT':'B',
     'ACGT':'N'
     }

BASES = {
    'A': 1, 'a': 1,
    'C': 2, 'c': 2,
    'G': 3, 'g': 3,
    'T': 4, 't': 4,
    'N': 0, 'n': 0,
    }

REVERSE_IUPAC = dict([(v,k) for (k,v) in IUPAC.items()])

try:
    from Bio.Seq import Seq
    from Bio.Alphabet import IUPAC as B_IUPAC
    from Bio import Translate
    alphabet = B_IUPAC.unambiguous_dna
    translator = Translate.unambiguous_dna_by_id[1]
except ImportError:
    pass

REVERSE_COMPL = {' ':'', 'A':'T', 'C':'G', 'G':'C', 'T':'A', '-':'-', 'N':'N', '[':']', ']':'[', 'a':'t', 'c':'g', 'g':'c', 't':'a', 'n':'n'}

REFSEQ_VERSION = 'hg18'

class CommandError(Exception):
    """ Raised when an external command-line program exits with an error code """
    pass


def timenow():
    """return current time as a string
    """
    return time.strftime('%d/%m/%Y %H:%M:%S', time.localtime(time.time()))


def reverse_compl(sequence):
    """ Reverse complement a genomic sequence
    """
    reverse = []

    ### Reverse the sequence, but keep any starting numeric values at
    ### the beginning.  I.e., the reverse-complement of the allele "7C"
    ### should be "7G", not "G7".
    if sequence[0].isdigit():
        for i, base in enumerate(sequence):
            if base.isdigit():
                reverse.append(base)
            else:
                bases = list(sequence[i:])
                break
    else:
        bases = list(sequence)
    bases.reverse()
    for base in bases:
        reverse.append(REVERSE_COMPL.get(base, base))
    return ''.join(reverse)

def strPercentage(numerator, denominator, format='%3.2f', ifZero='0.0'):
    """ Returns a string version of numerator/denominator, catching
        devide-by-zero errors and returning default
    """
    if denominator == 0:
        return ifZero
    else:
        percentage = (float(numerator)/float(denominator)) * 100.0
        return format % (percentage)
        
def getGenename():
    """ An algorithm to determine the name of the gene we
        are working on in this directory ...
    """
    ### If the file "HUGO.txt" exists, it should be considered
    ### canonical
    if os.path.exists('HUGO.txt'):
        hugofile = open('HUGO.txt', 'r')
        hugoname = hugofile.readline().strip()
        return hugoname
    else:
        ### Otherwise, the name of the current directory is
        ### the name of the gene
        top = os.getcwd()
        genename = os.path.basename(top)
        return genename

class LocusOrRange(object):
    """ Represents a chromosomal locus tuple (chrom, strand, posBefore, posAfter)
    """
    def __init__(self, lstring=None, chrom=None, strand=None, posBefore=None, posAfter=None):
        if lstring:
            self.parse(lstring)
        else:
            self.chrom = chrom
            self.strand = strand
            self.posBefore = posBefore
            self.posAfter = posAfter
            self.calcSize()
            self.calcAbspos()

    def __str__(self):
        size = self.size
        if not (self.chrom and self.strand and self.posBefore and self.posAfter):
            return ''
        elif size == 1:
            return '%s%s:%s' % (self.chrom, self.strand, self.posBefore+1)
        elif size > 1:
            return '%s%s:%s-%s' % (self.chrom, self.strand, self.posBefore+1, self.posAfter-1)
        elif size < 1:
            return '%s%s:%s^%s' % (self.chrom, self.strand, self.posBefore+1, self.posAfter-1)
        else:
            return '%s%s:%s-%s' % (self.chrom, self.strand, self.posBefore+1, self.posAfter-1)

    def calcSize(self):
        if self.posBefore and self.posAfter:
            self.size = self.posAfter - self.posBefore - 1
        else:
            self.size = None

    def calcAbspos(self):
        size = self.size
        if size == 1:
            self.abspos = int(self.posBefore+1)
        elif size > 1:
            self.abspos = '%s-%s' % (self.posBefore, self.posAfter)
        elif size < 1:
            self.abspos = '%s^%s' % (self.posBefore, self.posAfter)
        else:
            self.abspos = None
        
    def toStringNoStrand(self):
        """ Convert this locus to a string, without the strand designation
        """
        size = self.spanSize
        if size == 1:
            return '%s:%s' % (self.chrom, self.posBefore+1)
        elif size > 1:
            return '%s:%s-%s' % (self.chrom, self.posBefore, self.posAfter)
        else:
            return '%s:%s^%s' % (self.chrom, self.posBefore, self.posAfter)
            
    def parse(self, locus):
        """ Parse locus information from a string desgination
        """
        chrom, abspos = locus.split(':')
        chrom = chrom.strip()
        if chrom.endswith('+') or chrom.endswith('-'):
            self.chrom = chrom[:-1]
            self.strand = chrom[-1]
        else:
            self.chrom = chrom
            self.strand = ''

        if '-' in abspos:
            start, end = abspos.split('-')
            start = int(start)
            abspos = int(start)
            end = int(end)
            self.posBefore = start - 1
            self.posAfter = end + 1
        elif '^' in abspos:
            start, end = abspos.split('^')
            start = int(start)
            abspos = int(start)
            end = int(end)
            self.posBefore = start - 1
            self.posAfter = end + 1
        else:
            abspos = int(abspos)
            self.posBefore = abspos - 1
            self.posAfter = abspos + 1
        self.calcSize()
        self.calcAbspos()

class RefGene(object):
    """
    ### Here's what we are interested in, to the best of my understanding of
    ### gene structure:
    ###                                                         
    ###                                                               CDNA
    ###                            ____________________________________|___________________________________
    ###                           |                                                                        |
    ###                           |                                                                        |                                           
    ###                         -- EXON,NON-CODING                                              -- EXON,CODING
    ###                        |                                                               |
    ###                        |   -- EXON,CODING                                              |            -- EXON,NON-CODING
    ###                        |  |                                                            |           |
    ### PROMOTER     5'UTR     |  |              INTRON       EXON,CODING           INTRON     |           |          3'UTR
    ### <-------- _______________________________________________________________________________________________________________________
    ###          |             |  |           |           |                   |                |           |   |                         |
    ###          V             V  V           V           V                   V                V           V   V                         V
    ###          |              __|___________             ___________________                  ___________|___                          |
    ###          |             |  |           |           |                   |                |           |   |                         |
    ###          |             |  |           |           |                   |                |           |   |                         |
    ###          |-----------------------------------------------------------------------------------------------------------------------|
    ###          |             |  |           |           |                   |                |           |   |                         |
    ###          |             |__|___________|           |___________________|                |___________|___|                         |
    ###          |                |                                                                        |                             |
    ###   
    ###          txStart          cdsStart                                                                 cdsEnd                        txEnd
    ###
    ### For a particular position within a gene, we need to determine the functional
    ### status as one of: UTR, EXON/NON_CODING, EXON/CODING, INTRON.  There may also
    ### be some researchers who are interested in markers that are intronic, but are
    ### "close" to the Intron/Exon Boundary, but this is by no means a clear definition.
    ###
    ### N.B. The order of exons in exonStarts, exonEnds, and exonFrames depends on
    ### the strand.  Any functionality added here which depends on the order of
    ### exons will need to take this into account.  For now, the only thing that depends
    ### on the strand is whether a given position is in the 3' or 5' direction
    ###
    """
    def __init__(self, row):
        """ Instantiate using a row returned from a MySQLdb.cursor.DictCursor after querying from
            the hg18.refGene table
        """
        self.hugo = row['name2']
        self.isoform = row['name']
        self.strand = row['strand']
        self.chrom = row['chrom']
        self.txStart = int(row['txStart'])+1
        self.txEnd = int(row['txEnd'])
        self.cdsStart = int(row['cdsStart'])+1
        self.cdsEnd = int(row['cdsEnd'])
        self.exonStarts = [int(start)+1 for start in row['exonStarts'].tostring().rstrip(',').split(',')]
        self.exonEnds = [int(end) for end in row['exonEnds'].tostring().rstrip(',').split(',')]
        self.exons = zip(self.exonStarts, self.exonEnds)

        self.locus = '%s%s:%s-%s' % (self.chrom, self.strand, self.txStart, self.txEnd)
        self.abspos = '%s-%s' % (self.txStart, self.txEnd)

        ### Calculate the weird, old-ish "features" table
        self.features = []
        ### Store the coding segments of each exon as features
        for exonStart, exonEnd in self.exons:
            ### Forget exons that are entirely before cdsStart
            ### or after cdsEnd
            if exonStart <= self.cdsStart and exonEnd <= self.cdsStart:
                continue
            elif exonStart >= self.cdsEnd and exonEnd >= self.cdsEnd:
                continue
            ### Otherwise, add a CDS feature for the part of the
            ### exon contained in the range cdsStart-cdsEnd
            fStart = max(exonStart, self.cdsStart)
            fEnd = min(exonEnd, self.cdsEnd)
            self.features.append(('CDS', fStart+1, fEnd))
        self.features.append(('RNA', self.txStart+1, self.txEnd)) ### The "feature" of the entire RNA from txStart to txEnd
        self.features.append(('UTR', self.txStart+1, self.cdsStart)) ### There's the 5' UTR from txStart to cdsStart        
        self.features.append(('UTR', self.cdsEnd+1, self.txEnd)) ### Then the 3' UTR from cdsEnd to txEnd

    def __str__(self):
        return str(self.hugo)

    def getRelpos(self, abspos):
        """ Given an absolute genome position, determine its relative position
        """
        if '-' in str(abspos):
            posBefore, posAfter = abspos.split('-')
        elif '^' in str(abspos):
            posBefore, posAfter = abspos.split('^')
        elif str(abspos).isdigit():
            posBefore = int(abspos) - 1
            posAfter = int(abspos) + 1
        else:
            raise 'BadAbspos', abspos
        posBefore = int(posBefore)
        posAfter = int(posAfter)
        
        if self.strand == '-':
            relBefore = self.cdsEnd - posAfter
            relAfter =  self.cdsEnd - posBefore
        elif self.strand == '+':
            relBefore = posBefore - self.cdsStart
            relAfter = posAfter - self.cdsStart
        else:
            raise 'UnknownStrand', self.strand

        ##raise 'TEST', str((self.strand, posBefore, posAfter, self.cdsStart, self.cdsEnd, relBefore, relAfter))
    
        relposSize = abs(relAfter - relBefore) - 1
        if relposSize == 1:
            return relBefore + 1
        elif relposSize < 1:
            return '%s^%s' % (relBefore, relAfter)
        else:
            return '%s-%s' % (relBefore, relAfter)

    def overlapsCodingRange(self, abspos):
        """ Determine if the given abspos overlaps a coding sequence in this gene
        """ 
        if '-' in str(abspos):
            start, end = abspos.split('-')
            for pos in range(int(start), int(end)+1):
                if self.cdsStart <= int(pos) <= self.cdsEnd:
                    return True
            return False
        elif '^' in str(abspos):
            start, end = abspos.split('^')
            for pos in range(int(start), int(end)+1):
                if self.cdsStart <= int(pos) <= self.cdsEnd:
                    return True
            return False
        elif (self.cdsStart <= int(abspos) <= self.cdsEnd):
            return True
        return False
        
    def overlapsExon(self, abspos):
        """ Determine if given "abspos" overlaps an exon in this gene
        """
        for exonStart, exonEnd in self.exons:
            if '-' in str(abspos):
                start, end = abspos.split('-')
                for pos in range(int(start), int(end)+1):
                    if exonStart <= int(pos) <= exonEnd:
                        return True
                return False
            elif '^' in str(abspos):
                start, end = abspos.split('^')
                for pos in range(int(start), int(end)+1):
                    if self.cdsStart <= int(pos) <= self.cdsEnd:
                        return True
                return False
            elif exonStart <= int(abspos) <= exonEnd:
                return True
        return False
     
    def getRoles(self, abspos):
        roles = []

        isInCodingRange = self.overlapsCodingRange(abspos)
        isExonic = self.overlapsExon(abspos)

        if abspos < self.txStart:
            if abspos >= self.txStart-1000:
                roles.append('Promoter')
            else:
                roles.append('Upstream')
        elif abspos > self.txEnd:
            roles.append('Downstream')
        elif abspos < self.cdsStart:
            if isExonic:
                roles.append('Exon/Non-coding')
            else:
                if self.strand == '+':
                    roles.append('5\'UTR')
                else:
                    roles.append('3\'UTR')
        elif abspos > self.cdsEnd:
            if isExonic:
                roles.append('Exon/Non-coding')
            else:
                if self.strand == '+':                            
                    roles.append('3\'UTR')
                else:
                    roles.append('5\'UTR')
        elif isInCodingRange:
            if isExonic:
                roles.append('Exon/Coding')
            else:
                roles.append('Intron')
        else:
            raise 'UnknownRole', 'abspos "%s" does not fit in gene (tx=%s-%s, cds=%s-%s, exons=%s)' % (txStart, txEnd, cdsStart, cdsEnd, str(exons))

        return roles

    def getGeneticChanges(self, abspos, variants):
        """ Given the possible variants at chromosomal position abspos, determine
            any "interesting" genetic changes induced by each variant.  That would
            be any changes modifying a splice-site or interefering with the promoter
            region
        """
        import getGenomeSpan

        ### Stitch together the MRNA Sequence ...
        sequence = []
        changes = []
        relpos = self.getRelpos(abspos)
        sequence = getGenomeSpan.getGenomeSpan(self.chrom, self.txStart, self.txEnd)

        sequence2 = list(sequence) ### This makes a copy

        ### For each variant specified, determine if there is an amino
        ### acid change in the protein when changing from the reference
        ### base at that position, to this "mutant" base
        return []
        if relpos < 0:
            return 'Promoter?'
        refBase = sequence[int(relpos)]
        for variant in variants:
            ### Don't bother calculating a matching base, or deletions
            if variant in (refBase, '-'):
                continue
            sequence2[positionInSequence] = variant
            resp.append('NUCLEOTIDE=%s/%s' % (sequence[positionInSequence], sequence2[positionInSequence]))
            sequence = ''.join(sequence)
            sequence2 = ''.join(sequence2)
        return changes

    def getProteinChanges(self, abspos, variants):
        """
        """
        import getGenomeSpan

        ### Stitch together the MRNA Sequence ...
        sequence = []
        variantPositions = []        
        for exonStart, exonEnd in self.exons:
            ### Forget exons that are entirely before cdsStart
            ### or after cdsEnd
            if exonStart <= self.cdsStart and exonEnd <= self.cdsStart:
                continue
            elif exonStart >= self.cdsEnd and exonEnd >= self.cdsEnd:
                continue
            ### Otherwise, add a CDS feature for the part of the
            ### exon contained in the range cdsStart-cdsEnd
            fStart = max(exonStart, self.cdsStart)
            fEnd = min(exonEnd, self.cdsEnd)
            exon = getGenomeSpan.getGenomeSpan(self.chrom, fStart, fEnd)
            sequence.append(exon)

            ### Determine the relative position within the mrna sequence
            if fStart <= abspos <= fEnd:
                variantPositions.append(len(sequence) - (fEnd - abspos))

        changes = []
        warnings = Set()
        for positionInSequence in variantPositions:
            sequence = list(''.join(sequence))

            ### For each variant specified, determine if there is an amino
            ### acid change in the protein when changing from the reference
            ### base at that position, to this "mutant" base
            refBase = sequence[positionInSequence]
            for variant in variants:
                sequence2 = list(sequence)
                change = {}
                change['REF_BASE'] = refBase

                ### Don't bother calculating a matching base, or deletions
                if variant in (refBase, '-'):
                    continue
                sequence2[positionInSequence] = variant
                change['VARIANT_BASE'] = sequence2[positionInSequence]
                sequence = ''.join(sequence)
                sequence2 = ''.join(sequence2)

                mrna = Seq(sequence, alphabet)
                mrna2 = Seq(sequence2, alphabet)
                if self.strand == '-':
                    mrna = mrna.reverse_complement()
                    mrna2 = mrna2.reverse_complement()

                protein = translator.translate(mrna)
                protein2 = translator.translate(mrna2)

                aminopos = positionInSequence/3
                change['AMINOPOS'] = aminopos+1
                
                ### These two warning flags are here temporarily so we can get an idea
                ### of how often this happens.  We've noticed sequences of stitched-together
                ### coding sequences that are not divisible by three, so when they are
                ### translated into proteins, the remaining bases get truncated.  This is
                ### a general problem in that I don't know how best to handle that, and
                ### whether it's supposed to happen at all.  It is a more specific problem
                ### when the SNP of interest is actually one of the truncated bases!
                if len(mrna) % 3 != 0:
                    warnings.add('MRNA_LENGTH_NOT_3_MULTIPLE')
                if len(protein) <= aminopos:
                    warnings.add('PROTEIN_LENGTH_LESS_THAN_AMINOPOS')
                    change['PROTEIN_LENGTH'] = len(protein)
                    continue

                refAmino = protein[aminopos]
                varAmino = protein2[aminopos]
                if varAmino != refAmino:
                    ### Some of these are rather verbose and are only useful for debugging
                    #change['HEAD'] = ''.join([str(i)[-2:].ljust(3) for i in range(1, len(protein)+1)])
                    #change['PRO1'] = '  '.join(protein.tostring())
                    #change['PRO2'] = '  '.join(protein2.tostring())
                    #change['SEQ1'] = mrna.tostring()
                    #change['SEQ2'] = mrna2.tostring()
                    change['VARIANT_POS'] = positionInSequence
                    change['CODING_LENGTH'] = len(sequence)
                    change['PROTEIN_LENGTH'] = len(protein)
                    change['REF_AMINO'] = refAmino
                    change['VAR_AMINO'] = varAmino
                changes.append(change)
                    
        return changes, warnings

    
INDETERMINATE = None
TOP = 0
BOT = 1
def determineTopBot(sequence):
    """ Implements the Illumina TOP/BOT algorithm

        Returns a tuple:

        (TOP_OR_BOT, N, ALLELE_DESIGNATIONS)
    """
    
    bopen = sequence.index('[')
    bclose = sequence.index(']')
    flank5 = sequence[:bopen]
    flank3 = sequence[bclose+1:]
    variant = sequence[bopen+1:bclose]
    variants = sorted(list(Set([v.upper() for v in variant.split('/')])))

    maxn = min((len(flank5), len(flank3)))+1
    
    topOrBot = INDETERMINATE
    for n in range(maxn):
        if n == 0:
            if len(variants) != 2 or '-' in variants or max([len(v) for v in variants]) > 1:
                continue
            else:
                if 'A' in variants and 'T' not in variants:
                    topOrBot = TOP
                    break
                elif 'T' in variants and 'A' not in variants:
                    topOrBot = BOT
                    break
                else:
                    continue
        else:
            left = flank5[-n].upper()
            right = flank3[n-1].upper()

            if left in ('A', 'T') and right in ('C', 'G'):
                topOrBot = TOP
                break
            elif left in ('C', 'G') and right in ('A', 'T'):
                topOrBot = BOT
                break
            else:
                continue

    if topOrBot == INDETERMINATE:
        n = None

    ### Now determine allele designations as A/B
    ### This is only applicable if there are exactly two unique alleles
    ### *and* the alleles are both simple nucleotide
    if len(variants) != 2 or '-' in variants or max([len(v) for v in variants]) > 1:
        alleleDesignations = dict([(v, None) for v in variants])
    else:
        vseq = ''.join(variants)
        if vseq == 'AC':
            alleleDesignations = {'A':'A', 'C':'B'}
        elif vseq == 'AG':
            alleleDesignations = {'A':'A', 'G':'B'}
        elif vseq == 'AT':
            if topOrBot == TOP:
                alleleDesignations = {'A':'A', 'T':'B'}
            elif topOrBot == BOT:
                alleleDesignations = {'T':'A', 'A':'B'}
            else:
                alleleDesignations = {'A':None, 'T':None}
        elif vseq == 'CG':
            if topOrBot == TOP:
                alleleDesignations = {'C':'A', 'G':'B'}
            elif topOrBot == BOT:
                alleleDesignations = {'G':'A', 'C':'B'}
            else:
                alleleDesignations = {'C':None, 'G':None}
        elif vseq == 'CT':
            alleleDesignations = {'T':'A', 'C':'B'}
        elif vseq == 'GT':
            alleleDesignations = {'T':'A', 'G':'B'}
        else:
            alleleDesignations = dict([(v, None) for v in variants])

    return topOrBot, n, alleleDesignations

def lineCount(infile):
    """ As effificiently as possible, count the number of lines in a file """
    return int(commands.getoutput('/usr/bin/wc -l %s' % (infile)).split()[0])

MAX_IN_SUBQUERY_SIZE = 1000
def batchList(myList, batch_size=None, attr=None):
    """ Break myList into a tuple of tuples, each of size
        no greater than chunk_size.

        This is used when constructing a WHERE clause with
        an IN on a list > 1000 items.

        If attr is given, it means that the items in myList
        are dicts (or behave as such), and item[attr] should
        be in the resulting tuples instead of item itself.
    """
    ### As a default, we will just return the data as a single chunk
    if not batch_size:
        batch_size = len(myList)
    nGroups = len(myList)/batch_size
    if len(myList) % batch_size != 0:
        nGroups += 1
    for i in range(nGroups):
        if attr:
            chunk = tuple([item[attr] for item in myList[i::nGroups]])
        else:
            chunk = tuple(list(myList)[i::nGroups])
        yield chunk

def chunkList(myList, batch_size=MAX_IN_SUBQUERY_SIZE, attr=None):
    return list(batchList(myList, batch_size, attr))                

def runCommand(command, verbose=False, raiseOnError=True):
    """ A generic command-line runner with error handling
    """
    stat, out = commands.getstatusoutput(command)
    if verbose:
        print 'CWD: %s' % (os.getcwd())
        print 'CMD: %s' % (command)
        print 'OUT: %s' % (out)
    if int(stat) != 0:
        if raiseOnError:
            raise CommandError('\nCommand "%s" failed with code "%s".\nCWD=%s\nOutput was: \n"""\n%s\n"""' % (command % (locals()), stat, os.getcwd(), out))
        else:
            print 'CommandError: Command failed with code "%s", CWD=%s' % (stat, os.getcwd())

