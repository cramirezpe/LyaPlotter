# %% 
import os
import numpy as np 
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table, join
import glob
import logging
from lyacolore import utils
import logging
from functools import cached_property
from pathlib import Path
log = logging.getLogger(__name__)

class QuickQuasarsSim():
    sim_class = 'QuickQuasars'
    def __init__(self,id_, path,nside=16, pixel=0, name=None, file_type = 'gaussian', compression =True):
        self.id_ = id_
        self.__name__ = name if name else path[path.rfind('/')+1:]
        self.sim_path = path
        self.file_type = file_type
        self.nside = nside
        self.pixel = pixel
        self.compression = compression
        self.dirname = utils.get_dir_name(self.sim_path+'/spectra-{}'.format(self.nside), self.pixel)

    def __str__(self):
        return "{} sim. Id: {}\tName: {}\tPath: {}".format(self.sim_class,self.id_,self.__name__,self.sim_path)

    @classmethod
    def search(cls,path):
        sim_paths = set()
        for path in glob.glob(path + '/**/spectra-*-*.fits',recursive=True):
            sim_paths.add( str(Path(path).parents[3]) )
        return sorted(list(sim_paths), key=lambda folder: os.stat(folder).st_ctime)

    def get_redshift_files(self):
        truth_file = utils.get_file_name(self.dirname,'truth',self.nside,self.pixel,self.compression)
        zbest_file = utils.get_file_name(self.dirname, 'zbest', self.nside, self.pixel,self.compression)

        self.Truth = TruthFile(truth_file,self.file_type)
        self.Best = BestFile(zbest_file, self.file_type)
        return self.Truth, self.Best      

    def get_spectra(self,arm,lr_max=1200.):
        spectra_file = utils.get_file_name(self.dirname,'spectra',self.nside,self.pixel,self.compression)
        
        self.Spectra = Spectra(arm,spectra_file, self.Truth.z_qso, self.Best.z_qso, lr_max, self.file_type)
        return self.Spectra

class CoLoReSim():
    sim_class = 'CoLoRe'

    def __init__(self,id_,path,name=None, file_type='gaussian'):
        self.id_ = id_
        
        # Lets give it a name, if it is not provided I'll put the name of the folder
        if name:
            self.__name__ = name
        else:
            self.__name__ = path[path.rfind('/')+1:]

        self.sim_path = path
        try:
            self.param_file = sorted(glob.glob(self.sim_path + '/*.cfg'), key=os.path.getmtime)[0]
        except:
            self.param_file = None
        self.file_type = file_type    
        self.CoLoReFile = {}

    @staticmethod
    def search(path):
        sim_paths = []
        for path in glob.glob(path+'/**/out_srcs_s1_0.fits',recursive=True):
            sim_paths.append( path[:path.rfind('/')] )
        return sorted(sim_paths, key=lambda folder: os.stat(folder).st_ctime)

    def __str__(self):
        return "{} sim. Id: {}\tName: {}\tPath: {}".format(self.sim_class,self.id_,self.__name__,self.sim_path)

    def get_Sources(self,ifile=0, lr_max=1200.):
        file = self.sim_path + '/out_srcs_s1_{}.fits'.format(ifile)
        self.CoLoReFile[ifile] = CoLoReFile(file,lr_max,self.file_type)
        return self.CoLoReFile[ifile]

class LyaCoLoReSim():
    sim_class = 'LyaCoLoRe'
    def __init__(self, id_, path,nside =16, pixel=0, name=None, file_type='gaussian', compression=True):
        self.id_ = id_
        self.__name__ = name if name else path[path.rfind('/')+1:]
        self.sim_path = path
        self.file_type = file_type
        self.nside = nside
        self.pixel = pixel
        self.compression = compression
        self.dirname = utils.get_dir_name(self.sim_path,self.pixel)
        self.picca_files = {}

    def __str__(self):
        return "{} sim. Id: {}\tName: {}\tPath: {}".format(self.sim_class,self.id_,self.__name__,self.sim_path)

    @classmethod
    def search(cls,path):
        sim_paths = set()
        for path in glob.glob(path + '/**/master.fits',recursive=True):
            sim_paths.add( str(Path(path).parents[0]) )
        return sorted(list(sim_paths), key=lambda folder: os.stat(folder).st_ctime)

    def get_Transmission(self,lr_max=1200.):
        transmission_file = utils.get_file_name(self.dirname,'transmission',self.nside,self.pixel,self.compression)
        self.Transmission = Transmission(transmission_file,lr_max,self.file_type)
        return self.Transmission

    def get_GaussianCoLoRe(self,lr_max=1200.):
        gaussian_colore_file = utils.get_file_name(self.dirname,'gaussian-colore',self.nside,self.pixel,self.compression)
        self.GaussianCoLoRe = GaussianCoLoRe(gaussian_colore_file, lr_max, self.file_type)
        return self.GaussianCoLoRe

    def get_PiccaStyleFile(self,file_name,lr_max=1200):
        file = utils.get_file_name(self.dirname,file_name, self.nside, self.pixel, self.compression)
        self.picca_files[file_name] = PiccaStyleFile(file, lr_max,self.file_type,file_name)
        return self.picca_files[file_name]

class FileSkeleton:
    def __init__(self, file_path, file_type):
        self.file_path  = file_path
        self.hdulist    = fits.open(self.file_path)
        self.file_type  = file_type
        try:
            self.RA     = self.hdulist[1].data['RA']
            self.DEC    = self.hdulist[1].data['DEC']
        except:
            self.RA     = None
            self.DEC    = None 

    def plot_locations(self, ax=None, **kwargs):
        '''Plot the locations of the QSO, they should be incorporated in self.RA and self.DEC'''
        if not ax: fig, ax = plt.subplots()
      
        Plotter.plot_locations(ax,self.RA,self.DEC,**kwargs)
        return

    def plot_qso_dist(self, ax=None, **kwargs):
        '''Return the n(z) QSO distribution assuming that it is stored in the variable self.z_qso'''
        if not ax: fig, ax = plt.subplots()
        bins = np.linspace(1.0,4.0,50)

        Plotter.plot_dist(ax,self.z_qso, bins)
        ax.set_xlabel(r'$z$')
        ax.set_ylabel('# QSOs')
        return

    def plot_pdf(self, values, values_name=' ',ax=None):
        ''' Return the pdf of the values given. It makes uses of the variables self.z that must be defined'''
        if not ax: fig, ax = plt.subplots()

        # Plot the pdf of deltas in redshift bins.
        z_bins = [(0,1),(1,2),(2,3),(3,)]
        if self.file_type == 'gaussian':
            d_bins = np.linspace(-5,5,100)
        elif self.file_type == '2lpt':
            d_bins = np.linspace(-1,4,100)

        for i,zbin in enumerate(z_bins):
            if len(zbin)==2:
                w = (self.z > zbin[0]) * (self.z<zbin[1])
                label = r'${}<z<{}$'.format(zbin[0], zbin[1])
            else:
                w = ((self.z>zbin[0]))
                label = r'${}<z$'.format(zbin[0])
            Plotter.plot_dist(ax,values=np.ravel(values[:,w]),bins=d_bins,weights=np.ravel(self.mask[:,w]),density=True,label=label)

        ax.set_xlabel('$\\delta$')
        ax.set_ylabel('$P(\delta)$')
        ax.legend()
        return

    def single_skewer(self, values_array, axis_values=None, value_name='', ax=None, mockid=None, **kwargs):
        '''Plot a single skewer from the FileClass. It returns the axis values and the values plotted. 
        
        Arguments:
        values_array -- All the values (for all skewers) for the value we want to compute. 

        Keyword arguments:
        axis_values -- Values we want to use as the x axis.
        value_name  -- For labelling purposes
        ax          -- Make the plot in an existing axis
        mockid      -- Id of the mock to plot (It should be defined correctly in self.id)   
        '''
        if not np.any(axis_values): axis_values = self.z
        if not ax:                  fig, ax = plt.subplots()
        
        if not mockid:              ind = 0
        else:                       ind = list(self.id).index(mockid)
        
        values = values_array[ind,:]

        return Plotter.plot_skewer(ax, axis_values=axis_values, values=values, value_name=value_name ,weights=self.mask[ind],**kwargs) 
         

    def mean_all_skewers(self, values_array, axis_values=None, value_name='', ax=None, kwargs_hline=None, **kwargs):
        '''Plot the mean value over all the skewers given. It applies the mask given by self.mask. It returns the axis values and the values plotted. 

        Arguments:
        values_array -- Values of all the skewers that we want to compute. 

        Keyword arguments:
        axis_values -- Values we want to use as the x axis.
        value_name  -- For labelling purposes
        ax          -- Make the plot in an existing axis
        kwargs_hline-- Style for horizonal line
        **kwargs    -- Style for the skewers line
        '''
        if not np.any(axis_values): axis_values = self.z
        if not ax:                  fig, ax = plt.subplots()

        return Plotter.plot_mean_all_skewers(ax,axis_values,values_array,value_name,self.mask, kwargs_hline=kwargs_hline, **kwargs)
        

    def std_all_skewers(self, values_array, axis_values=None, value_name='', ax=None, kwargs_hline=None, **kwargs):
        '''Plot the std over all the skewers given. It applies the mask given by self.mask. It returns the axis values and the values plotted. 

        Arguments:
        values_array -- Values of all the skewers that we want to compute. 

        Keyword arguments:
        axis_values -- Values we want to use as the x axis.
        value_name  -- For labelling purposes
        ax          -- Make the plot in an existing axis
        kwargs_hline-- Style for horizonal line
        **kwargs    -- Style for the skewers line
        '''
        if not np.any(axis_values): axis_values = self.z
        if not ax:                  fig, ax = plt.subplots()

        return Plotter.plot_std_all_skewers(ax,axis_values, values_array, value_name, self.mask, kwargs_hline=kwargs_hline, **kwargs)
        

class CoLoReFile(FileSkeleton):
    def __init__(self, file_path, lr_max, file_type):
        FileSkeleton.__init__(self,file_path,file_type)

        self.z_qso = self.hdulist[1].data['Z_COSMO']
        self.N_qso = len(self.z_qso)

        self.z = np.asarray( self.hdulist[4].data['Z'] )
        self.wavelength = utils.lya_rest * (1+self.z)

        self.mask = utils.make_IVAR_rows(lr_max,self.z_qso,np.log10(utils.lya_rest*(1+self.z))) 
        self.delta_skewers = self.hdulist[2].data
        self.vrad = self.hdulist[3].data
        self.lr_max = lr_max
        self.id = list(range(self.N_qso))

class TruthFile(FileSkeleton):
    def __init__(self,file_path, file_type):
        FileSkeleton.__init__(self,file_path, file_type)
        self.z_qso = self.hdulist[1].data['Z']

class BestFile(FileSkeleton):
    def __init__(self,file_path,file_type):
        FileSkeleton.__init__(self,file_path, file_type)
        self.z_qso = self.hdulist[1].data['Z']

class Spectra(FileSkeleton):
    def __init__(self,arm,spectra_file,z_qso_truth,z_qso_best,lr_max,file_type):
        assert arm.lower() in ['r','b','z']
        self.arm    = arm.lower()
        FileSkeleton.__init__(self,spectra_file,file_type)
        self.lr_max = lr_max

        self.z_qso_best     = z_qso_best
        self.z_qso_truth    = z_qso_truth

        #self.z_qso = self.z_qso_best
        self.z_qso  = self.z_qso_truth

        self.RA     = self.hdulist[1].data['TARGET_RA']
        self.DEC    = self.hdulist[1].data['TARGET_DEC']

        if self.arm == 'r':
            self.flux = self.hdulist['R_FLUX'].data
            self.wavelength = self.hdulist['R_WAVELENGTH'].data
        elif self.arm == 'b':
            self.flux = self.hdulist['B_FLUX'].data
            self.wavelength = self.hdulist['B_WAVELENGTH'].data
        else:
            self.flux = self.hdulist['Z_FLUX'].data
            self.wavelength = self.hdulist['Z_WAVELENGTH'].data

        self.z = (self.wavelength - utils.lya_rest)/utils.lya_rest
        self.mask = utils.make_IVAR_rows( lr_max, self.z_qso, np.log10(self.wavelength) )

        self.id = self.hdulist[1].data['TARGETID']

class Transmission(FileSkeleton):
    def __init__(self, file_path, lr_max, file_type):
        FileSkeleton.__init__(self,file_path,file_type)
        self.lr_max = lr_max

        self.z_qso = self.hdulist[1].data['Z']
        self.z_qso_noRSD = self.hdulist[1].data['Z_noRSD']
        self.N_qso = len(self.z_qso)
        

        self.wavelength = self.hdulist[2].data
        self.z = (self.wavelength - utils.lya_rest)/utils.lya_rest

        # Exctract Lyman absorption
        self.lya_absorption = self.hdulist['F_LYA'].data
        self.delta_lya_absorption = self.lya_absorption/Computations.mean_per_pixel(self.lya_absorption,mask=False)

        self.lyb_absorption = self.hdulist['F_LYB'].data

        # Setting mask
        self.mask = utils.make_IVAR_rows(lr_max,self.z_qso,np.log10(utils.lya_rest*(1+self.z)))

        self.id = self.hdulist[1].data['MOCKID']

class GaussianCoLoRe(FileSkeleton):
    def __init__(self,file_path,lr_max=1200., file_type = 'gaussian'):
        FileSkeleton.__init__(self, file_path, file_type)

        self.z_qso = self.hdulist[1].data['Z_COSMO']
        self.N_qso = len(self.z_qso)

        self.z = np.asarray( self.hdulist[4].data['Z'])
        self.wavelength = utils.lya_rest * (1+self.z)

        self.mask = utils.make_IVAR_rows(lr_max,self.z_qso,np.log10(utils.lya_rest*(1+self.z)))
        self.delta_skewers = self.hdulist[2].data
        self.vrad = self.hdulist[3].data
        
        self.id = self.hdulist[1].data['MOCKID']
        self.lr_max = lr_max

class PiccaStyleFile(FileSkeleton):
    def __init__(self,file_path,lr_max=1200., file_type = 'gaussian', name = ''):
        FileSkeleton.__init__(self,file_path,file_type)
        self.__name__ = name 

        self.z_qso = self.hdulist[3].data['Z']
        self.N_qso = len(self.z_qso)

        self.wavelength = 10**self.hdulist[2].data
        self.z = (self.wavelength - utils.lya_rest)/utils.lya_rest

        self.mask = utils.make_IVAR_rows(lr_max,self.z_qso,np.log10(utils.lya_rest*(1+self.z)))
        self.values = self.hdulist[0].data.transpose()

        self.id = self.hdulist[3].data['THING_ID']
        self.lr_max = lr_max

class Plotter:
    @staticmethod
    def plot_locations(ax,RA,DEC,**kwargs):
        phi = RA*np.pi/180
        theta = np.pi/2 - DEC*np.pi/180
        ax.scatter(phi/np.pi,np.cos(theta),**kwargs)
        ax.set_xlim(0.0,2.0)
        ax.set_ylim(-1.0,1.0)
        ax.set_xlabel(r'$\phi/\pi$')
        ax.set_ylabel(r'$\cos(\theta)$')
        return phi, theta
    
    @staticmethod
    def plot_dist(ax,values,bins,histtype= 'step', **kwargs):
        ax.hist(values,bins,histtype =histtype, **kwargs)
        return bins,values

    @staticmethod
    def plot_skewer(ax, axis_values,values,value_name,weights=None,print_info=False, **kwargs):
        if np.any(weights):
            mean_values = np.average(values,weights=weights)
            mean_values_squared = np.average(values**2,weights=weights)
            sigma_values = np.sqrt(mean_values_squared-mean_values**2)

            if print_info:
                print('Mean {} = {}'.format(value_name,mean_values))
                print('Std {} = {}'.format(value_name, sigma_values))

        ax.plot(axis_values, values, **kwargs)
        return axis_values, values
    
    @staticmethod
    def plot_mean_all_skewers(ax,axis_values,values,value_name,mask,print_info=True, kwargs_hline = None, **kwargs):
        overall_mean = Computations.overall_mean(values,mask)
        if print_info: print('Mean over all pixels = {}'.format(overall_mean))
        kwargs_hline = kwargs_hline or {'c': 'grey'}
        
        w = mask.sum(axis=0)>0
        axis_values = axis_values[w]
        values = Computations.mean_per_pixel(values,mask)

        ax.plot(axis_values, values,**kwargs)
        ax.axhline(y=overall_mean, **kwargs_hline)
        return axis_values, values 
        
    @staticmethod
    def plot_std_all_skewers(ax,axis_values,values,value_name,mask,print_info=True, kwargs_hline = None, **kwargs):
        overall_sigma = Computations.overall_sigma(values,mask)
        if print_info: print('Std over all pixels = {:1.4f}'.format(overall_sigma))
        kwargs_hline = kwargs_hline or {'c': 'grey'}

        w           = mask.sum(axis=0)>0
        axis_values = axis_values[w]
        values      = Computations.std_per_pixel(values,mask)


        ax.plot(axis_values, values, **kwargs)
        ax.axhline(y=overall_sigma, **kwargs_hline)
        return axis_values, values

class Computations:
    @classmethod
    def overall_mean(cls,values,mask):
        return np.average(values, weights=mask)
    
    @classmethod
    def overall_sigma(cls,values,mask):
        overall_mean_squared = np.average(values**2, weights=mask)
        return np.sqrt(overall_mean_squared - cls.overall_mean(values,mask)**2)

    @classmethod
    def mean_per_pixel(cls,values,mask=False):
        if not np.any(mask):
            return np.average(values,axis=0)
        else:
            w = mask.sum(axis=0)>0 # This is because one value is not masked I should put the axis value
            return np.average(values[:,w],weights=mask[:,w],axis=0)
    
    @classmethod
    def std_per_pixel(cls,values,mask):
        w = mask.sum(axis=0)>0
        mean = cls.mean_per_pixel(values,mask)
        mean_squared = cls.mean_per_pixel(values**2,mask)
        return np.sqrt(mean_squared-mean**2)