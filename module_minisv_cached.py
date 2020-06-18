from module_lyacolore_multifile_cached import (cached_property, FilesSkeleton)
import numpy as np
from astropy.io import fits

class zBestFilesMiniSV(FilesSkeleton):
    def __init__(self, file_paths):
        FilesSkeleton.__init__(self,file_paths=file_paths, parent_sim=None)

    @cached_property
    def qsos_target_ids(self):
        tids = []
        for h in self.hdulists:
            # h_tids = h[1].data['TARGETID']
            # mask = (h[1].data['SPECTYPE'][:].astype(str)=='QSO')&(h[1].data['ZWARN'][:]==0)
            # tids.append(h_tids[mask])
            h_tids = h[2].data['TARGETID']
            mask = h[2].data['CMX_TARGET'][:].astype(str)=='4096'
            tids.append(h_tids[mask])
        return tids

    @cached_property
    def RA(self):
        return  np.concatenate( [hdulist[2].data['TARGET_RA'][idx]  for hdulist,idx in zip(self.hdulists,self.idx_qsos_2)] )

    @cached_property
    def DEC(self):
        return  np.concatenate( [hdulist[2].data['TARGET_DEC'][idx] for hdulist,idx in zip(self.hdulists,self.idx_qsos_2)] )

    @cached_property
    def z_qso(self):
        return  np.concatenate( [hdulist[1].data['Z'][idx] for hdulist,idx in zip(self.hdulists, self.idx_qsos_1)] )

    @cached_property
    def N_qso(self):
        return len(self.z_qso)

    @cached_property
    def flux_r(self):
        return  np.concatenate( [hdulist[2].data['FLUX_R'][idx]     for hdulist,idx in zip(self.hdulists,self.idx_qsos_2)] )

    @cached_property
    def flux_g(self):
        return  np.concatenate( [hdulist[2].data['FLUX_G'][idx]     for hdulist,idx in zip(self.hdulists,self.idx_qsos_2)] )

    @cached_property
    def flux_z(self):
        return  np.concatenate( [hdulist[2].data['FLUX_Z'][idx]     for hdulist,idx in zip(self.hdulists,self.idx_qsos_2)] )

    @cached_property
    def idx_qsos_1(self):
        idx_qsos_1 = []
        for hdulist, qsos_ids in zip( self.hdulists, self.qsos_target_ids ):
            tid_1 = hdulist[1].data['TARGETID']
            idx_qsos_1.append( self.get_qso_indx_from_target_ids(tid_1, qsos_ids) )
        return idx_qsos_1

    @cached_property
    def idx_qsos_2(self):
        idx_qsos_2 = []
        for hdulist, qsos_ids in zip( self.hdulists, self.qsos_target_ids ):
            tid_2 = hdulist[2].data['TARGETID']
            idx_qsos_2.append( self.get_qso_indx_from_target_ids(tid_2, qsos_ids) )
        return idx_qsos_2
        
    @staticmethod
    def get_qso_indx_from_target_ids(thids_dest, qsos_ids):
        try:
            indxs = []
            for tid in qsos_ids:
                indxs.append( np.where(thids_dest==tid )[0][0] )
            return indxs
        except IndexError as e:
            print(thids_source)
            raise


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

