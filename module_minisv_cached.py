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
            # h_tids_1 = h[1].data['TARGETID']
            # mask_1 = (h[1].data['SPECTYPE'][:].astype(str)=='QSO')&(h[1].data['ZWARN'][:]==0)
            # tids_1 = set(h_tids_1[mask_1])

            # h_tids_2 = h[2].data['TARGETID']
            # mask_2 = h[2].data['CMX_TARGET'][:].astype(str)=='4096'
            # tids_2 = set(h_tids_2[mask_2])

            # tids.append(list(tids_2.intersection(tids_1) ) )

            # -------------------------------------

            h_tids = h[2].data['TARGETID']
            mask = h[2].data['CMX_TARGET'][:].astype(str)=='4096'
            tids.append(h_tids[mask])

            # --------------------------------------

            # h_tids  = h[1].data['TARGETID']
            # mask    = (h[1].data['SPECTYPE'][:].astype(str)=='QSO')&(h[1].data['ZWARN'][:]==0)  
            # tids.append(h_tids[masks])
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
        except IndexError:
            print(qsos_ids)
            raise

class DrqFile(FilesSkeleton):
    def __init__(self, file_paths):
        FilesSkeleton.__init__(self, file_paths=file_paths)

    @cached_property
    def z_qso(self):
        return  np.concatenate( [hdulist[1].data['Z']       for hdulist in self.hdulists] )

    @cached_property
    def N_qso(self):
        return len(self.z_qso)

    def fiberid(self):
        return np.concatenate( [hdulist[1].data['fiberid']   for hdulist in self.hdulists] )

    @cached_property
    def thingid(self):
        return np.concatenate( [hdulist[1].data['thing_id']   for hdulist in self.hdulists] )

    @cached_property
    def plate(self):
        return np.concatenate( [hdulist[1].data['plate']     for hdulist in self.hdulists] )
    
    @cached_property
    def best_night(self):
        return np.concatenate( [hdulist[1].data['mjd']       for hdulist in self.hdulists] )

    @cached_property
    def coadd_object(self,andes_location='/global/cfs/cdirs/desi/spectro/redux/andes/tiles/'):
        paths = set()
        for plate,night in zip(map(str,self.plate),map(str,self.best_night)):
            tile = plate[:-1]
            petal    = plate[-1]
            paths.add(f'{andes_location}{tile}/{night}/coadd-{petal}-{tile}-{night}.fits')
        return CoaddFilesMiniSV(paths)

    @cached_property
    def match_coadd_files(self):
        idx_coadd = []
        idx_drq = []
        y = self.coadd_object

        for i,(plate, night, tid) in enumerate( zip(map(str,self.plate), self.best_night,self.thingid) ) :
            tile = plate[:-1]
            petal = int( plate[-1] )
            logic = np.logical_and.reduce((y.targetid==tid,y.night == night,y.tile_spec==tile,y.petal_spec ==petal))
            try:
                match = np.where(logic)[0][0]
                idx_drq.append(i)
                idx_coadd.append(match)
            except:
                print('Error: qso in catalog not in data')
                raise

        return idx_drq, idx_coadd

    @property
    def idx_drq(self):
        return self.match_coadd_files[0]
        
    @property
    def idx_coadd(self):
        return self.match_coadd_files[1]

    @cached_property
    def flux_r(self):
        return  self.coadd_object.flux_r[self.idx_coadd]

    @cached_property
    def flux_g(self):
        return  self.coadd_object.flux_g[self.idx_coadd]

    @cached_property
    def flux_z(self):
        return  self.coadd_object.flux_z[self.idx_coadd]



class CoaddFilesMiniSV(FilesSkeleton):
    def __init__(self, file_paths):
        FilesSkeleton.__init__(self,file_paths=file_paths)
    
    @cached_property
    def RA(self):
        return np.concatenate( [hdulist[1].data['TARGET_RA']    for hdulist in self.hdulists] )
    
    @cached_property
    def DEC(self):
        return np.concatenate( [hdulist[1].data['TARGET_DEC']   for hdulist in self.hdulists] ) 

    @cached_property
    def targetid(self):
        return np.concatenate( [hdulist[1].data['TARGETID']     for hdulist in self.hdulists] ) 

    @property
    def tile_spec(self):
        return np.concatenate( [np.full(hdulist[1].header['NAXIS2'],path.split('-')[-2]) for hdulist,path in zip(self.hdulists,self.file_paths)] )
        # return np.asarray( [path.split('-')[-2] for path in self.file_paths] )
        
    @property
    def petal_spec(self):
        return np.concatenate( [hdulist[1].data['PETAL_LOC']    for hdulist in self.hdulists] )

    @property
    def night(self):
        return np.concatenate( [hdulist[1].data['NIGHT'] for hdulist in self.hdulists] )

    @property
    def flux_r(self):
        return  np.concatenate( [hdulist[1].data['FLUX_R']   for hdulist in self.hdulists] )

    @property
    def flux_g(self):
        return  np.concatenate( [hdulist[1].data['FLUX_G']   for hdulist in self.hdulists] )

    @property
    def flux_z(self):
        return  np.concatenate( [hdulist[1].data['FLUX_Z']   for hdulist in self.hdulists] )