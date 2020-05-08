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
from pathlib import Path
log = logging.getLogger(__name__)

class cached_property(object):
    """
    Descriptor (non-data) for building an attribute on-demand on first use.
    """
    def __init__(self, factory):
        """
        <factory> is called such: factory(instance) to build the attribute.
        """
        self._attr_name = factory.__name__
        self._factory = factory

    def __get__(self, instance, owner):
        # Build the attribute.
        attr = self._factory(instance)

        # Cache the value; hide ourselves.
        setattr(instance, self._attr_name, attr)

        return attr

class QuickQuasarsSim():
    sim_class = 'QuickQuasars'

    def __init__(self,id_, path, nside=16, name=None, file_type = 'gaussian', compression =True):
        self.id_        = id_
        self.__name__   = name if name else path[path.rfind('/')+1:]
        self.sim_path   = path
        self.file_type  = file_type
        self.nside      = nside
        self.compression= compression

    def __str__(self):
        return "{} sim. Id: {}\tName: {}\tPath: {}".format(self.sim_class,self.id_,self.__name__,self.sim_path)

    @classmethod
    def search(cls,path):
        sim_paths = set()
        for path in glob.glob(path + '/**/spectra-*-*.fits',recursive=True):
            sim_paths.add( str(Path(path).parents[3]) )
        return sorted(list(sim_paths), key=lambda folder: os.stat(folder).st_ctime)

    def get_spectra(self, arm, pixels=[0], lr_max=1200., redshift_to_use= 'truth'):
        ''' Get spectra for a single arm:

        Arguments:
        arm             -- Arm for which we want to get spectra (r,b or z)
        pixels          -- Pixels we want to use
        lr_max          -- Set maximum wavelength for mask
        redhisft_to_use -- We need a redshift for the Quasar to compute quantities (truth or best)
        '''
        check_is_list(pixels)
        truth_files  = []
        zbest_files  = []
        spectra_files = []
        for pixel in pixels:
            dirname     = utils.get_dir_name(self.sim_path+'/spectra-{}'.format(self.nside), pixel )
            truth_files.append(  utils.get_file_name(dirname,'truth',  self.nside,pixel,self.compression) )
            zbest_files.append(  utils.get_file_name(dirname,'zbest',  self.nside,pixel,self.compression) )
            spectra_files.append(utils.get_file_name(dirname,'spectra',  self.nside,pixel,self.compression) )

        return Spectra(arm, spectra_files, truth_files, zbest_files, pixels, lr_max, redshift_to_use, self)      

class CoLoReSim():
    sim_class = 'CoLoRe'

    def __init__(self,id_,path,name=None, file_type='gaussian'):
        self.id_ = id_
        
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

    @staticmethod
    def search(path):
        sim_paths = []
        for path in glob.glob(path+'/**/out_srcs_s1_0.fits',recursive=True):
            sim_paths.append( path[:path.rfind('/')] )
        return sorted(sim_paths, key=lambda folder: os.stat(folder).st_ctime)

    def __str__(self):
        return "{} sim. Id: {}\tName: {}\tPath: {}".format(self.sim_class,self.id_,self.__name__,self.sim_path)

    def get_Sources(self, ifiles=[0], lr_max=1200.):
        ''' Get sources from a CoLoRe simulations:

        Arguments
        ifiles          -- Array with the srcs files that we want to use. 
        lr_max          -- Set maximum wavelength for mask
        '''
        check_is_list(ifiles)
        files            = [self.sim_path + '/out_srcs_s1_{}.fits'.format(ifile) for ifile in ifiles]

        return CoLoReFiles(files, lr_max, self)
        # self.CoLoReFiles = CoLoReFiles(files,lr_max,self.file_type)
        # return self.CoLoReFiles

class LyaCoLoReSim():
    sim_class = 'LyaCoLoRe'
    def __init__(self, id_, path, nside=16, name=None, file_type='gaussian', compression=True):
        self.id_        = id_
        self.__name__   = name if name else path[path.rfind('/')+1:]
        self.sim_path   = path
        self.file_type  = file_type
        self.nside      = nside
        self.compression= compression
        self.picca_files= {}

    def __str__(self):
        return "{} sim. Id: {}\tName: {}\tPath: {}".format(self.sim_class,self.id_,self.__name__,self.sim_path)

    @classmethod
    def search(cls,path):
        sim_paths = set()
        for path in glob.glob(path + '/**/master.fits',recursive=True):
            sim_paths.add( str(Path(path).parents[0]) )
        return sorted(list(sim_paths), key=lambda folder: os.stat(folder).st_ctime)

    def get_Transmission(self,  pixels=[0], lr_max=1200.):
        ''' Get transmission for the given pixels:

        Arguments:
        pixels          -- Pixels we want to use
        lr_max          -- Set maximum wavelength for mask
        '''
        check_is_list(pixels)
        files = []
        for pixel in pixels:
            dirname = utils.get_dir_name(self.sim_path, pixel)
            files.append( utils.get_file_name(dirname, 'transmission',     self.nside, pixel, self.compression) )

        return Transmission(files, lr_max, self)
        # self.Transmission = Transmission(files, lr_max, self.file_type)
        # return self.Transmission

    def get_GaussianCoLoRe(self,pixels=[0], lr_max=1200.):
        ''' Get colore skewers/velocity for the given pixels:

        Arguments:
        pixels          -- Pixels we want to use
        lr_max          -- Set maximum wavelength for mask
        '''
        check_is_list(pixels)
        files = []
        for pixel in pixels:
            dirname = utils.get_dir_name(self.sim_path, pixel)
            files.append(  utils.get_file_name(dirname, 'gaussian-colore', self.nside, pixel, self.compression)   )
        return GaussianCoLoRe(files, lr_max, self)
        # self.GaussianCoLoRe = GaussianCoLoRe(files, lr_max, self.file_type)
        # return self.GaussianCoLoRe

    def get_PiccaStyleFiles(self,file_name, pixels=[0], lr_max=1200):
        ''' Get values stored in file file_name. It should be the format 'picca' files that we can get from LyaCoLoRe:

        Arguments:
        file_name       -- Name of the file (e.g. picca-flux-noRSD-notnorm)
        pixels          -- Pixels we want to use
        lr_max          -- Set maximum wavelength for mask
        '''
        check_is_list(pixels)
        files = []
        for pixel in pixels:
            dirname = utils.get_dir_name(self.sim_path, pixel)
            files.append( utils.get_file_name(dirname, file_name,          self.nside, pixel, self.compression) )
        return PiccaStyleFiles(files, lr_max, file_name, self)
        # self.picca_files[file_name] = PiccaStyleFiles(files, lr_max,self.file_type,file_name)
        # return self.picca_files[file_name]

class FilesSkeleton:
    def __init__(self, file_paths, parent_sim):
        check_is_list(file_paths)

        self.sim        = parent_sim

        self.file_paths = file_paths
        self.hdulists   = [fits.open(path) for path in self.file_paths]
        
    @cached_property
    def RA(self):
        return  np.concatenate( [hdulist[1].data['RA'] for hdulist in self.hdulists] )

    @cached_property
    def DEC(self):
        return  np.concatenate( [hdulist[1].data['DEC'] for hdulist in self.hdulists] )

    def plot_locations(self, ax=None, **kwargs):
        '''Plot the locations of the QSO, they should be incorporated in self.RA and self.DEC
        
        Arguments:
        ax              -- Set the axis where to plot it. By default it will create a new axis. 
        kwargs          -- All unmatched kwargs will be sent to the plotter.
        '''
        if not ax: fig, ax = plt.subplots()
      
        Plotter.plot_locations(ax,self.RA,self.DEC,**kwargs)
        return

    def plot_qso_dist(self, ax=None, **kwargs):
        '''Return the n(z) QSO distribution assuming that it is stored in the variable self.z_qso
        
        Arguments:
        ax              -- Set the axis where to plot it. By default it will create a new axis. 
        kwargs          -- All unmatched kwargs will be sent to the plotter.
        '''
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
        if self.sim.file_type == 'gaussian':
            d_bins = np.linspace(-5,5,100)
        elif self.sim.file_type == '2lpt':
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
        kwargs      -- All unmatched kwargs will be sent to the plotter.
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
        

class CoLoReFiles(FilesSkeleton):
    def __init__(self, file_paths, lr_max, parent_sim):
        FilesSkeleton.__init__(self,file_paths,parent_sim)        
        self.lr_max             = lr_max

    @cached_property
    def z_qso(self):
        return  np.concatenate( [hdulist[1].data['Z_COSMO'] for hdulist in self.hdulists] )

    @cached_property
    def N_qso(self):
        return  len(self.z_qso)
    
    @cached_property
    def id(self):
        return  list(range(self.N_qso)) # Apparently they do not have an Id so I'll let their name just be a number

    @cached_property
    def z(self):
        return  np.asarray( self.hdulists[0][4].data['Z'])
    
    @cached_property
    def wavelength(self):
        return  utils.lya_rest * (1 + self.z)

    @cached_property
    def mask(self):
        return  utils.make_IVAR_rows(self.lr_max, self.z_qso, np.log10(utils.lya_rest*(1+self.z)))

    @cached_property
    def delta_skewers(self):
        return  np.vstack( [hdulist[2].data for hdulist in self.hdulists] )
            
    @cached_property
    def vrad(self):
        return  np.vstack( [hdulist[3].data for hdulist in self.hdulists] )

class TruthFiles(FilesSkeleton):
    def __init__(self, file_paths, parent_sim):
        FilesSkeleton.__init__(self, file_paths, parent_sim)

class BestFiles(FilesSkeleton):
    def __init__(self, file_paths, parent_sim):
        FilesSkeleton.__init__(self, file_paths, parent_sim)

class Spectra(FilesSkeleton):
    def __init__(self, arm, spectra_files, truth_files, zbest_files, pixels, lr_max, redshift_to_use, parent_sim):
        FilesSkeleton.__init__(self,spectra_files, parent_sim)
        self.truth      = TruthFiles(truth_files,  parent_sim)
        self.zbest      = BestFiles(zbest_files,   parent_sim)

        assert arm.lower() in ['r','b','z']
        assert redshift_to_use.lower() in ['truth','best']
        self.arm             = arm.lower()
        self.lr_max          = lr_max
        self.redshift_to_use = redshift_to_use.lower()
       
    @cached_property
    def z_qso_best(self):
        return np.concatenate( [hdulist[1].data['Z'] for hdulist in self.zbest.hdulists] )

    @cached_property
    def z_qso_truth(self):
        return np.concatenate( [hdulist[1].data['Z'] for hdulist in self.truth.hdulists] )

    @cached_property
    def z_qso(self):
        if self.redshift_to_use == 'best':
            return self.z_qso_best
        else:
            return self.z_qso_truth

    @cached_property
    def N_qso(self):
        return
    
    @cached_property
    def id(self):
        return  np.concatenate( [hdulist[1].data['TARGETID']   for hdulist in self.hdulists] )

    @cached_property
    def z(self):
        return  (self.wavelength - utils.lya_rest)/utils.lya_rest
    
    @cached_property
    def wavelength(self):
        return  self.hdulists[0]['{}_WAVELENGTH'.format(self.arm.upper())].data

    @cached_property
    def mask(self):
        return  utils.make_IVAR_rows( self.lr_max, self.z_qso, np.log10(self.wavelength) )

    @cached_property
    def RA(self):
        return  np.concatenate( [hdulist[1].data['TARGET_RA']  for hdulist in self.hdulists] )

    def DEC(self):
        return  np.concatenate( [hdulist[1].data['TARGET_DEC'] for hdulist in self.hdulists] )

    @cached_property
    def flux(self):
        return  np.vstack( [hdulist['{}_FLUX'.format(self.arm.upper())].data for hdulist in self.hdulists ] )

class Transmission(FilesSkeleton):
    def __init__(self, file_paths, lr_max, parent_sim):
        FilesSkeleton.__init__(self,file_paths, parent_sim)
        self.lr_max         = lr_max

    @cached_property
    def z_qso(self):
        return  np.concatenate( [hdulist[1].data['Z']        for hdulist in self.hdulists] )

    @cached_property
    def z_qso_noRSD(self):
        return  np.concatenate( [hdulist[1].data['Z_noRSD']  for hdulist in self.hdulists] )

    @cached_property
    def N_qso(self):
        return  len(self.z_qso)
    
    @cached_property
    def id(self):
        return  np.concatenate( [hdulist[1].data['MOCKID']   for hdulist in self.hdulists] )

    @cached_property
    def z(self):
        return  (self.wavelength - utils.lya_rest)/utils.lya_rest
    
    @cached_property
    def wavelength(self):
        return  self.hdulists[0][2].data

    @cached_property
    def mask(self):
        return  utils.make_IVAR_rows(self.lr_max,self.z_qso,np.log10(self.wavelength))

    @cached_property
    def lya_absorption(self):
        return  np.vstack(  [hdulist['F_LYA'].data for hdulist in self.hdulists] )

    @cached_property
    def lyb_absorption(self):
        return  np.vstack(  [hdulist['F_LYB'].data for hdulist in self.hdulists] )

    @cached_property
    def delta_lya_absorption(self):
        return  self.lya_absorption/Computations.mean_per_pixel(self.lya_absorption,mask=False)
    
    @cached_property
    def delta_lyb_absorption(self):
        return  self.lyb_absorption/Computations.mean_per_pixel(self.lyb_absorption,mask=False)

class GaussianCoLoRe(FilesSkeleton):
    def __init__(self,file_paths, lr_max, parent_sim):
        FilesSkeleton.__init__(self, file_paths, parent_sim)
        self.lr_max             = lr_max

    @cached_property
    def z_qso(self):
        return  np.concatenate( [ hdulist[1].data['Z_COSMO']  for hdulist in self.hdulists] )

    @cached_property
    def N_qso(self):
        return  len(self.z_qso)
    
    @cached_property
    def id(self):
        return  np.concatenate( [ hdulist[1].data['MOCKID']   for hdulist in self.hdulists] )

    @cached_property
    def z(self):
        return  self.hdulists[0][4].data['Z']
    
    @cached_property
    def wavelength(self):
        return   utils.lya_rest * (1+self.z)

    @cached_property
    def mask(self):
        return  utils.make_IVAR_rows(self.lr_max,self.z_qso,np.log10(utils.lya_rest*(1+self.z)))

    @cached_property
    def delta_skewers(self):
        return  np.vstack( [hdulist[2].data for hdulist in self.hdulists] )

    @cached_property
    def vrad(self):
        return  np.vstack( [hdulist[3].data for hdulist in self.hdulists] )
        
class PiccaStyleFiles(FilesSkeleton):
    def __init__(self,file_paths,lr_max, name, parent_sim):
        FilesSkeleton.__init__(self,file_paths, parent_sim)
        self.__name__       = name 
        self.lr_max         = lr_max

    @cached_property
    def z_qso(self):    
        return  np.concatenate( [ hdulist[3].data['Z']        for hdulist in self.hdulists] )

    @cached_property
    def N_qso(self):
        return  len(self.z_qso)
    
    @cached_property
    def id(self):
        return  np.concatenate( [ hdulist[3].data['THING_ID'] for hdulist in self.hdulists] )

    @cached_property
    def z(self):
        return  (self.wavelength - utils.lya_rest)/utils.lya_rest
    
    @cached_property
    def wavelength(self):
        return  10**self.hdulists[0][2].data

    @cached_property
    def mask(self):
        return  utils.make_IVAR_rows(self.lr_max,self.z_qso,np.log10(self.wavelength))

    @cached_property
    def values(self):
        tmp = np.vstack( [hdulist[0].data   for hdulist in self.hdulists] )
        return tmp.transpose()
        
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

def check_is_list(x):
    if isinstance(x,list) or isinstance(x,tuple):
        return True
    else:
        type_ = type(x)
        raise ValueError('A list/tuple was expected. Type of {} was {}'.format(x,type_),x)