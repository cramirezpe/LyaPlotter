from module_lyacolore_multifile_cached import (cached_property, FilesSkeleton)
import numpy as np
from astropy.io import fits

class zBestFilesMiniSV(FilesSkeleton):
    def __init__(self, file_paths):
        FilesSkeleton.__init__(self,file_paths=file_paths, parent_sim=None)

    @cached_property
    def select_qsos(self):
        return [(h[1].data['SPECTYPE'][:].astype(str)=='QSO')&(h[1].data['ZWARN'][:]==0) for h in self.hdulists ]

    @cached_property
    def idx_qsos_2(self):
        # In the hdulist 2 there are duplicated values, hence instead of a select we should apply an index
        idx_qsos_2 = []
        tids1 = self.tid
        tids2 = self.tid2

        for i in range(len(self.hdulists)):
            idx_qsos_2.append( [] )
            for tid in tids1[i]:
                idx_qsos_2[i].append( np.where(tids2[i] == tid)[0][0] )
        
        return idx_qsos_2

    @cached_property
    def tid(self):
        return [hdulist[1].data['TARGETID'][select] for hdulist,select in zip(self.hdulists, self.select_qsos) ]

    @cached_property
    def tid2(self):
        return [ hdulist[2].data['TARGETID'] for hdulist in self.hdulists ]
        
    @cached_property
    def RA(self):
        return  np.concatenate( [hdulist[2].data['TARGET_RA'][idx] for hdulist,idx in zip(self.hdulists,self.idx_qsos_2)] )

    @cached_property
    def DEC(self):
        return  np.concatenate( [hdulist[2].data['TARGET_DEC'][idx] for hdulist,idx in zip(self.hdulists,self.idx_qsos_2)] )

    # @cached_property
    # def RA(self):
    #     return np.concatenate( [hdulist[2].data['TARGET_RA']  for hdulist in self.hdulists] )
    
    # @cached_property
    # def DEC(self):
    #     return np.concatenate( [hdulist[2].data['TARGET_DEC'] for hdulist in self.hdulists] ) 


class DrqFile(FilesSkeleton):
    def __init__(self, file_paths):
        FilesSkeleton.__init__(self, file_paths=file_paths)


class CoaddFilesMiniSV(FilesSkeleton):
    def __init__(self, file_paths):
        FilesSkeleton.__init__(self,file_paths=file_paths)
    
    @cached_property
    def RA(self):
        return np.concatenate( [hdulist[1].data['TARGET_RA']  for hdulist in self.hdulists] )
    
    @cached_property
    def DEC(self):
        return np.concatenate( [hdulist[1].data['TARGET_DEC'] for hdulist in self.hdulists] ) 

