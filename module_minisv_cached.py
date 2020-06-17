from module_lyacolore_multifile_cached import (cached_property, FilesSkeleton)
import numpy as np
from astropy.io import fits

class zBestFilesMiniSV(FilesSkeleton):
    def __init__(self, file_paths):
        FilesSkeleton.__init__(self,file_paths=file_paths, parent_sim=None)


    @cached_property
    def RA(self):
        return np.concatenate( [hdulist[2].data['TARGET_RA']  for hdulist in self.hdulists] )
    
    @cached_property
    def DEC(self):
        return np.concatenate( [hdulist[2].data['TARGET_DEC'] for hdulist in self.hdulists] ) 

    @cached_property
    def z_qso(self):
        return np.concatenate( [hdulist[1].data['Z']          for hdulist in self.hdulists] )
    
    @cached_property
    def N_qso(self):
        return len(self.z_qso)

class DrqFile(FilesSkeleton):
    def __init__(self, file_paths):
        FilesSkeleton.__init__(self, file_paths=file_paths)

    @cached_property
    def z_qso(self):
        return  np.concatenate( [hdulist[1].data['Z'] for hdulist in self.hdulists] )

    @cached_property
    def N_qso(self):
        return len(self.z_qso)

class CoaddFilesMiniSV(FilesSkeleton):
    def __init__(self, file_paths):
        FilesSkeleton.__init__(self,file_paths=file_paths)
    
    @cached_property
    def RA(self):
        return np.concatenate( [hdulist[1].data['TARGET_RA']  for hdulist in self.hdulists] )
    
    @cached_property
    def DEC(self):
        return np.concatenate( [hdulist[1].data['TARGET_DEC'] for hdulist in self.hdulists] ) 

