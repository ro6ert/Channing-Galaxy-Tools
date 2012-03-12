import os
import shutil
from subprocess import Popen

import pysam

def make_destination_folder(destination):
    out = os.path.abspath(destination)
    out_dir = os.path.dirname(out)
    if not os.path.exists(out_dir): os.makedirs(out_dir)


def SamBamFactory(filename):
    if filename.endswith('.sam'): return Sam(filename)
    elif filename.endswith('.bam'): return Bam(filename)
    else:
        raise ValueError, 'Can not tell which file this is.'


class SamBamFile(object):
    def __init__(self, filename):
        self.filename = os.path.abspath(filename)



class Sam(SamBamFile):
    @classmethod
    def file_extension(cls):
        return 'sam'

    def to_bam(self, destination):
        make_destination_folder(destination)
        out = os.path.abspath(destination)

        # Run samtools
        process = Popen('samtools view -h -b -S -o %s %s' % (out, self.filename), shell=True)
        return_value = process.wait()

        # Things running well, we return the created file
        if not return_value: return Bam(out)

    def to_sam(self, destination):
        out_full_path = os.path.abspath(destination)
        if not (self.filename == out_full_path):
            shutil.copy(self.filename, out_full_path)

        return Sam(out_full_path)


    def toggled_format(self, destination):
        return self.to_bam(destination)
    
    def same_format(self, destination):
        return self.to_sam(destination)
            

class Bam(SamBamFile):
    
    @classmethod
    def file_extension(cls):
        return 'bam'

    
    def to_bam(self, destination):
        out_full_path = os.path.abspath(destination)
        if not (self.filename == out_full_path):
            shutil.copy(self.filename, out_full_path)

        return Bam(out_full_path)

    def to_sam(self, destination):
        make_destination_folder(destination)
        out = os.path.abspath(destination)

        # Run samtools
        process = Popen('samtools view -h -o %s %s' % (out, self.filename), shell=True)
        return_value = process.wait()

        # Things running well, we return the created file
        if not return_value: return Sam(out)

    def toggled_format(self, destination):
        return self.to_sam(destination)

    def same_format(self, destination):
        return self.to_bam(destination)


        


if __name__ == '__main__':
    s = Sam('kwdfkjsdfg.sam')
    print s.file_extension()
