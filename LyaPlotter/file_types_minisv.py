'''
    Module build to analyse different data related to the Andes Data Release.
'''
from LyaPlotter.file_types import *
import numpy as np



class DrqFile(FilesBase):
    ''' class to handle data from catalogue files.

    It uses the object CoaddFilesMiniSV to obtain information only available in coadd files. It uses the hacked best_night inside the field 'mjd' to match the correct coadd file.
    '''

    def __init__(self, file_paths):
        FilesBase.__init__(self, file_paths=file_paths)

    @cached_property
    def z(self):
        with self.open_hdulists():
            return  np.concatenate( [hdulist[1].data['Z']       for hdulist in self.hdulists] )

    @cached_property
    def N_obj(self):
        return len(self.z)

    def fiberid(self):
        with self.open_hdulists():
            return np.concatenate( [hdulist[1].data['fiberid']   for hdulist in self.hdulists] )

    @cached_property
    def thingid(self):
        with self.open_hdulists():
            return np.concatenate( [hdulist[1].data['thing_id']   for hdulist in self.hdulists] )

    @cached_property
    def plate(self):
        with self.open_hdulists():
            return np.concatenate( [hdulist[1].data['plate']     for hdulist in self.hdulists] )
    
    @cached_property
    def best_night(self):
        with self.open_hdulists():
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



class CoaddFilesMiniSV(FilesBase):
    def __init__(self, file_paths):
        FilesBase.__init__(self,file_paths=file_paths)
    
    @cached_property
    def RA(self):
        with self.open_hdulists():
            return np.concatenate( [hdulist[1].data['TARGET_RA']    for hdulist in self.hdulists] )
    
    @cached_property
    def DEC(self):
        with self.open_hdulists():
            return np.concatenate( [hdulist[1].data['TARGET_DEC']   for hdulist in self.hdulists] ) 

    @cached_property
    def targetid(self):
        with self.open_hdulists():
            return np.concatenate( [hdulist[1].data['TARGETID']     for hdulist in self.hdulists] ) 

    @cached_property
    def tile_spec(self):
        with self.open_hdulists():
            return np.concatenate( [np.full(hdulist[1].header['NAXIS2'],path.split('-')[-2]) for hdulist,path in zip(self.hdulists,self.file_paths)] )
        # return np.asarray( [path.split('-')[-2] for path in self.file_paths] )
        
    @cached_property
    def petal_spec(self):
        with self.open_hdulists():
            return np.concatenate( [hdulist[1].data['PETAL_LOC']    for hdulist in self.hdulists] )

    @cached_property
    def night(self):
        with self.open_hdulists():
            return np.concatenate( [hdulist[1].data['NIGHT'] for hdulist in self.hdulists] )

    @property
    def flux_r(self):
        with self.open_hdulists():
            return  np.concatenate( [hdulist[1].data['FLUX_R']   for hdulist in self.hdulists] )

    @property
    def flux_g(self):
        with self.open_hdulists():
            return  np.concatenate( [hdulist[1].data['FLUX_G']   for hdulist in self.hdulists] )

    @property
    def flux_z(self):
        with self.open_hdulists():
            return  np.concatenate( [hdulist[1].data['FLUX_Z']   for hdulist in self.hdulists] )


class zBestFilesMiniSV(FilesBase):
    ''' Class to handle data from MiniSV zbest files.

    The main difference from other FilesBase is the fact that here a filter is applied to eliminate all the object that are not QSOs.

    It is needed to define correctly the method qsos_target_ids in order to select (for each fits file) the targetids that we will consider. 

    Attributes (extending FilesBase attributes):
        qsos_target_ids (array of arrays of ints): Containing the target ids corresponding to quasars for each fits file provided.
        idx_qsos_1 (array of arrays of ints): Indices corresponding to quasars for the first header for all the fits files.
        idx_qsos_2 (array of arrays of ints): Indices corresponding to quasars for the second header for ll the fits files. 

    '''
    def __init__(self, file_paths):
        FilesBase.__init__(self,file_paths=file_paths, parent_sim=None)

    @cached_property
    def qsos_target_ids(self):
        with self.open_hdulists():
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
        with self.open_hdulists():
            return  np.concatenate( [hdulist[2].data['TARGET_RA'][idx]  for hdulist,idx in zip(self.hdulists,self.idx_qsos_2)] )

    @cached_property
    def DEC(self):
        with self.open_hdulists():
            return  np.concatenate( [hdulist[2].data['TARGET_DEC'][idx] for hdulist,idx in zip(self.hdulists,self.idx_qsos_2)] )

    @cached_property
    def z(self):
        with self.open_hdulists():
            return  np.concatenate( [hdulist[1].data['Z'][idx] for hdulist,idx in zip(self.hdulists, self.idx_qsos_1)] )

    @cached_property
    def tile(self):
        with self.open_hdulists():
            return np.concatenate( [hdulist[2].data['TILEID'][idx] for hdulist,idx in zip(self.hdulists, self.idx_qsos_2)] )

    @cached_property
    def flux_r(self):
        with self.open_hdulists():
            return  np.concatenate( [hdulist[2].data['FLUX_R'][idx]     for hdulist,idx in zip(self.hdulists,self.idx_qsos_2)] )

    @cached_property
    def flux_g(self):
        with self.open_hdulists():
            return  np.concatenate( [hdulist[2].data['FLUX_G'][idx]     for hdulist,idx in zip(self.hdulists,self.idx_qsos_2)] )

    @cached_property
    def flux_z(self):
        with self.open_hdulists():
            return  np.concatenate( [hdulist[2].data['FLUX_Z'][idx]     for hdulist,idx in zip(self.hdulists,self.idx_qsos_2)] )

    @cached_property
    def idx_qsos_1(self):
        with self.open_hdulists():
            idx_qsos_1 = []
            for hdulist, qsos_ids in zip( self.hdulists, self.qsos_target_ids ):
                tid_1 = hdulist[1].data['TARGETID']
                idx_qsos_1.append( self.get_qso_indx_from_target_ids(tid_1, qsos_ids) )
            return idx_qsos_1

    @cached_property
    def idx_qsos_2(self):
        with self.open_hdulists():
            idx_qsos_2 = []
            for hdulist, qsos_ids in zip( self.hdulists, self.qsos_target_ids ):
                tid_2 = hdulist[2].data['TARGETID']
                idx_qsos_2.append( self.get_qso_indx_from_target_ids(tid_2, qsos_ids) )
            return idx_qsos_2

    @staticmethod
    def get_qso_indx_from_target_ids(thids_dest, qsos_ids):
        '''Obtain indices where the corresponding target_id is on qsos_ids.
        '''
        try:
            indxs = []
            for tid in qsos_ids:
                indxs.append( np.where(thids_dest==tid )[0][0] )
            return indxs
        except IndexError:
            print(qsos_ids)
            raise