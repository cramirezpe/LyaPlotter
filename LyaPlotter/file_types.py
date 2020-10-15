'''
    Module build to structure the different file types from different origins.

    The objective of this module is to homogenize the way data is obtained from different types of files. 

    data_dictionary is defined for each class it consists on a dictionary where it is stored the information on how to search the information inside the fits file. It is composed of four variables. 
        - name: Name of the card to read
        - field: Name of the field to read
        - vstack: Whether to vstack the data or not
        - only_one: Only grab one element (for shape issues)
'''

from lyacolore import utils
import matplotlib.pyplot as plt
import healpy as hp
import logging
from astropy.io import fits
import numpy as np
from LyaPlotter.plotter import Plotter
from LyaPlotter.computations import Computations
from contextlib import contextmanager


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

class FilesBase:
    '''Base class to handle output files

    Attributes:
        file_paths (list of str): List of paths to the different output files.
        sim (Sim object): Object of the sim the output files are coming from.
        N_files (int): Number of output files.
        hdulists (fits): Fits from the output files. 

        RA (list of float): RA from files.
        DEC (list of float): DEC from files.
    '''

    data_dictionary = {
        'RA'    : (1,   'RA',   False,  False),
        'DEC'   : (1,   'DEC',  False,  False),
        'z'     : (1,   'Z',    False,  False),
    }


    def __init__(self, file_paths, parent_sim=None, downsampling =1):
        '''Inits the File class setting the different information

        Args:
            file_paths (list of str): List of paths to the different output files.
            parent_sim (Sim object, optional): Object of the sim the output files are coming from.
            downsampling (float): Downsampling to apply to the data.

        '''
        if isinstance(file_paths, str): #pragma: no cover
            file_paths = [file_paths]

        self.sim        = parent_sim
        self.downsampling= downsampling

        self.file_paths = file_paths
        self.N_files    = len(file_paths)

    @contextmanager
    def open_hdulists(self):
        '''Context manager to open fits files and close it afterwards

        Args:
            file_paths (list of str): Paths to the different fits files.
            
        Returns:
            List of all the (opened) fits files to handle them.
        '''
        try:
            self.hdulists = [fits.open(path) for path in self.file_paths]
            yield self.hdulists
        except: #pragma: no cover
            print('Reading fits files failed:',self.file_paths)
            raise
        finally:
            [x.close() for x in self.hdulists]

    def get_data(self, name, field=None, vstack=False, only_one=False):
        with self.open_hdulists():
            if only_one:
                if not field:
                    return np.asarray( self.hdulists[0][name].data )
                else:
                    return np.asarray( self.hdulists[0][name].data[field])

            if self.downsampling != 1:
                output = list()
                for i,hdulist in enumerate(self.hdulists):
                    np.random.seed(i)
                    length = hdulist[name].header['NAXIS2']
                    if not field: 
                        new_data = hdulist[name].data[np.random.choice(a=[True,False], size=length, p = [self.downsampling,1-self.downsampling])]
                    else:
                        new_data = hdulist[name].data[field][np.random.choice(a=[True,False], size=length, p = [self.downsampling,1-self.downsampling])]                
                    output.append(new_data)            
            else:
                if not field:
                    output = [hdulist[name].data for hdulist in self.hdulists]
                else:
                    output = [hdulist[name].data[field] for hdulist in self.hdulists]   

            if vstack: 
                return np.vstack( output ) 
            else:
                return np.concatenate( output )   

    @cached_property
    def RA(self):
        '''Cached property obtaining RA from files

        Returns:
            A list of RA values (floats)
        '''
        return self.get_data((*self.data_dictionary['RA']))

    @cached_property
    def DEC(self):
        '''Cached property obtaining DEC from files

        Returns:
            A list of DEC values (floats)
        '''
        return self.get_data((*self.data_dictionary['DEC']))

    @cached_property
    def z(self):
        return self.get_data((*self.data_dictionary['z']))
            
   
    @cached_property
    def N_obj(self):
        return  len(self.z)

    def plot_locations(self, ax=None, **kwargs): #pragma: no cover
        '''Plot the locations from the file, they should be incorporated in self.RA and self.DEC
        
        Args:
            ax (axis object, optional): Set the axis where to plot it. By default it will create a new axis. 
            **kwargs (optional): Additional arguments to the plot.
        '''
        if not ax: fig, ax = plt.subplots()
        
        Plotter.plot_locations(ax,self.RA,self.DEC,**kwargs)
        return

    def plot_footprint(self, nside=64, ax=None, **kwargs): #pragma: no cover
        '''Plot the locations of the QSO in a heatmap, they should be incorporated in self.RA and self.DEC
        
        Args:
            nside (int, obptional): Pixelization of the plot.
            **kwargs (optional): Additional arguments to the plot.
        '''
        if not ax: 
            log.info('Axis not provided')
            fig, ax = plt.subplots()
        else:
            log.info('Axis provided')
        npix = hp.nside2npix(nside)
        pix  = hp.ang2pix(nside, self.RA, self.DEC,     lonlat=True)
        n = np.bincount(pix, minlength=npix)

        plt.axes(ax)
        hp.mollview(n, hold=True, unit='# objects', **kwargs)
        if not ax: plt.show()
        return

    def plot_dist(self, ax=None, bins=None, **kwargs): #pragma: no cover
        '''Return the n(z) distribution of objects assuming that it is stored in the variable self.z
        
        Arguments:
            ax (axis object, optional): Set the axis where to plot it. By default it will create a new axis. 
            **kwargs (): Additional arguments to the plot.
        '''
        if not ax: fig, ax = plt.subplots()
        if bins is None: bins = np.linspace(1.0,4.0,50)

        Plotter.plot_dist(self.z, bins,ax=ax, **kwargs)
        ax.set_xlabel(r'$z$')
        ax.set_ylabel('# QSOs')
        return

class FilesSkewerBase(FilesBase):
    '''Extension of the FilesBase class to allow skewer analysis.

    Attributes (extending FilesBase attributes):
        z_skewer (array of float): Redshift position of each skewer.
        wavelength (array of float): Wavelength of the LyAlpha line at the corresponding z_skewer redshift.
        mask (array of bool): Mask applied to the skewers.
        id (array of int): Id of each QSO skewer.

    '''
    wavelength = None
    mask       = None

    @cached_property
    def delta_skewers(self):
        return self.get_data(*self.data_dictionary['delta_skewers'])

    @cached_property
    def z_skewer(self):
        '''Obtain the pixelization of the fits file (in z)

        This is the same for every fits file so I should only grab it from one case
        '''
        return self.get_data(*self.data_dictionary['z_skewer'])

    @cached_property
    def wavelength(self):
        '''Obtain the pixelization of the fits file (in z)

        This is the same for every fits file so I should only grab it from one case
        '''
        return self.get_data(*self.data_dictionary['wavelength'])

    @cached_property
    def id(self):
        return self.get_data((*self.data_dictionary['id']))

    def single_skewer(self, values_array, axis_values=None, value_name='', ax=None, mockid=None, **kwargs): #pragma: no cover
        '''Plot a single skewer from the FileClass. It returns the axis values and the values plotted. 
        
        Args:
            values_array (2dim array): All the values (for all skewers) for the value we want to compute. 
            axis_values (array of float, optional): Values we want to use as the x axis. (By default self.z)
            value_name (str, optional): Name of the variable (for labeling)
            ax (axis object, optional): Set the axis where to plot it. By default it will create a new axis.
            mockid (int, optional): Id of the mock to plot (It should be defined correctly in self.id)
            **kwargs (optional): Additional arguments to the plot.

        Returns:
            The axis values and the values of the plot.
        '''
        if not np.any(axis_values): axis_values = self.z_skewer
        if not ax:                  fig, ax = plt.subplots()
        
        if not mockid:              ind = 0
        else:                       ind = list(self.id).index(mockid)
        
        values = values_array[ind,:]

        return Plotter.plot_skewer(ax, axis_values=axis_values, values=values, value_name=value_name ,weights=self.mask[ind],**kwargs) 
         

    def mean_all_skewers(self, values_array, axis_values=None, value_name='', ax=None, kwargs_hline=None, **kwargs): #pragma: no cover
        '''Plot the mean value over all the skewers given. It applies the mask given by self.mask. It returns the axis values and the values plotted. 

        Args:
            values_array (2dim array): All the values (for all skewers) for the value we want to compute. 
            axis_values (array of float, optional): Values we want to use as the x axis (By default self.z).
            value_name (str, optional): Name of the variable (for labeling)
            ax (axis object, optional): Set the axis where to plot it. By default it will create a new axis.
            **kwargs_hline (optional): Arguments to the hline plot. 
            **kwargs (optional): Additional arguments to the plot.

        Returns:
            The axis values and the values of the plot.
        '''
        if not np.any(axis_values): axis_values = self.z_skewer
        if not ax:                  fig, ax = plt.subplots()

        return Plotter.plot_mean_all_skewers(ax,axis_values,values_array,value_name,self.mask, kwargs_hline=kwargs_hline, **kwargs)
        

    def std_all_skewers(self, values_array, axis_values=None, value_name='', ax=None, kwargs_hline=None, **kwargs): #pragma: no cover
        '''Plot the std over all the skewers given. It applies the mask given by self.mask. It returns the axis values and the values plotted. 

        Args:
            values_array (2dim array): All the values (for all skewers) for the value we want to compute. 
            axis_values (array of float, optional): Values we want to use as the x axis (By default self.z).
            value_name (str, optional): Name of the variable (for labeling)
            ax (axis object, optional): Set the axis where to plot it. By default it will create a new axis.
            **kwargs_hline (optional): Arguments to the hline plot. 
            **kwargs (optional): Additional arguments to the plot.

        Returns:
            The axis values and the values of the plot.
        '''
        if not np.any(axis_values): axis_values = self.z_skewer
        if not ax:                  fig, ax = plt.subplots()

        return Plotter.plot_std_all_skewers(ax,axis_values, values_array, value_name, self.mask, kwargs_hline=kwargs_hline, **kwargs)

    def plot_pdf(self, values, values_name=' ',ax=None, **kwargs): #pragma: no cover
        ''' Return the pdf of the values given. It makes uses of the variables self.z_skewer that must be defined'''
        if not ax: fig, ax = plt.subplots()

        # Plot the pdf of deltas in redshift bins.
        z_bins = [(0,1),(1,2),(2,3),(3,)]
        if self.sim.file_type == 'gaussian':
            d_bins = np.linspace(-5,5,100)
        elif self.sim.file_type == '2lpt':
            d_bins = np.linspace(-1,4,100)

        for i,zbin in enumerate(z_bins):
            if len(zbin)==2:
                w = (self.z_skewer > zbin[0]) * (self.z_skewer<zbin[1])
                label = r'${}<z<{}$'.format(zbin[0], zbin[1])
            else:
                w = ((self.z_skewer>zbin[0]))
                label = r'${}<z$'.format(zbin[0])
            Plotter.plot_dist(ax=ax,values=np.ravel(values[:,w]),bins=d_bins,weights=np.ravel(self.mask[:,w]),density=True,label=label, **kwargs)

        ax.set_xlabel('$\\delta$')
        ax.set_ylabel('$P(\delta)$')
        ax.legend()
        return


class LyaCoLoReMasterFile(FilesBase):
    ''' Class to handle LyaCoLoRe master files.
    '''
    data_dictionary = {
    'RA'            : (1,   'RA',       False,  False),
    'DEC'           : (1,   'DEC',      False,  False),
    'z_noRSD'       : (1,   'Z_QSO_NO_RSD', False,  False),
    'z'             : (1,   'Z_QSO_RSD',False,  False),
    'id'            : (1,   'MOCKID',   False,  False),
    'pixnum'        : (1,   'PIXNUM',   False,  False),
    'filenum'       : (1,   'FILENUM',  False,  False)
    }

    def __init__(self, file_paths):
        FilesBase.__init__(self, file_paths=file_paths)

    @cached_property
    def id(self):
        return self.get_data((*self.data_dictionary['id']))
    
    @cached_property
    def z_noRSD(self):
        return self.get_data((*self.data_dictionary['z_noRSD']))

    @cached_property
    def pixnum(self):
        return self.get_data((*self.data_dictionary['pixnum']))

    @cached_property
    def filenum(self):
        return self.get_data((*self.data_dictionary['filenum']))

class CoLoReFiles(FilesSkewerBase):
    '''Class to handle CoLoRe output files. 

    Attributes (extending FilesSkewerBase):
        lr_max (float): Maximum wavelength for masking.
    '''

    data_dictionary = {
        'RA'            : (1,   'RA',       False,  False),
        'DEC'           : (1,   'DEC',      False,  False),
        'E1'            : (1,   'E1',       False,  False),
        'E2'            : (1,   'E2',       False,  False),
        'z'             : (1,   'Z_COSMO',  False,  False),
        'delta_skewers' : (2,   None,       True,   False),
        'vrad'          : (3,   None,       True,   False),
        'z_skewer'      : (4,   'Z',        False,  True)
    }

    def __init__(self, file_paths, lr_max, parent_sim,downsampling):
        '''Inits the CoLoReFiles object

        Args (extending FilesBase init args):
            lr_max (float): Maximum wavelength for masking.
        '''
        FilesSkewerBase.__init__(self,file_paths,parent_sim, downsampling)        
        self.lr_max             = lr_max

    @cached_property
    def id(self): #pragma: no cover
        '''Id of each skewer

        Apparently they do not have an Id so I'll let their name just be a number
        '''
        return  list(range(self.N_obj))

    @cached_property
    def E1(self):
        return self.get_data(*self.data_dictionary['E1'])
    
    @cached_property
    def E2(self):
        return self.get_data(*self.data_dictionary['E2'])

    @cached_property
    def wavelength(self):
        return  utils.lya_rest * (1 + self.z_skewer)

    @cached_property
    def mask(self):
        return  utils.make_IVAR_rows(self.lr_max, self.z, np.log10(utils.lya_rest*(1+self.z_skewer)))

    @cached_property
    def vrad(self):
        return self.get_data(*self.data_dictionary['vrad'])

class TruthFiles(FilesBase):
    '''TruthFiles can be handled with the information given in FilesBase'''
    pass

class BestFiles(FilesBase):
    '''BestFiles can be handled with the information given in FilesBase'''
    pass

class TransmissionFiles(FilesSkewerBase):
    '''Class to handle LyaCoLoRe transmission files.

    Attributes (extending FilesSkewerBase):
        lr_max (float): Maximum wavelength for masking.
    '''

    data_dictionary = {
        'RA'            : (1,       'RA',       False,  False),
        'DEC'           : (1,       'DEC',      False,  False),
        'z'             : (1,       'Z',        False,  False),
        'Z_noRSD'       : (1,       'Z_noRSD',  False,  False),
        'id'            : (1,       'MOCKID',   False,  False),
        'delta_skewers' : (2,       None,       True,   False),
        'vrad'          : (3,       None,       True,   False),
        'lya_absorption': ('F_LYA', None,       True,   False),
        'lyb_absorption': ('F_LYB', None,       True,   False),
        'wavelength'    : (2,       None,       True,   True)
    }

    def __init__(self, file_paths, lr_max, parent_sim, downsampling):
        '''Inits the TransmissionFiles object

        Args (extending FilesBase init args):
            lr_max (float): Maximum wavelength for masking.
        '''
        FilesBase.__init__(self,file_paths, parent_sim, downsampling)
        self.lr_max         = lr_max

    @cached_property
    def z_noRSD(self):
        return self.get_data(*self.data_dictionary['Z_noRSD'])
    
    @cached_property
    def z_skewer(self):
        return  (self.wavelength - utils.lya_rest)/utils.lya_rest
    
    @cached_property
    def mask(self):
        return  utils.make_IVAR_rows(self.lr_max,self.z,np.log10(self.wavelength))

    @cached_property
    def lya_absorption(self):  
        return self.get_data(*self.data_dictionary['lya_absorption'])
 
    @cached_property
    def lyb_absorption(self):
        return self.get_data(*self.data_dictionary['lyb_absorption'])

    @cached_property
    def delta_lya_absorption(self):
        return  self.lya_absorption/Computations.mean_per_pixel(self.lya_absorption,mask=False)
    
    @cached_property
    def delta_lyb_absorption(self):
        return  self.lyb_absorption/Computations.mean_per_pixel(self.lyb_absorption,mask=False)

    @cached_property
    def vrad(self):
        return self.get_data(*self.data_dictionary['vrad'])

class GaussianCoLoReFiles(CoLoReFiles):
    '''Class to handle GaussianCoLoRe files output by LyaCoLoRe. 

    It behaves exactly as the CoLoReFiles class but the data dictionary should be changed. (Also the id can be recovered from the fits file now) 
    '''
    data_dictionary = {
        'RA'            : (1,   'RA',       False,  False),
        'DEC'           : (1,   'DEC',      False,  False),
        'z'             : (1,   'Z_COSMO',  False,  False),
        'id'            : (1,   'MOCKID',   False,  False),
        'delta_skewers' : (2,   None,       True,   False),
        'vrad'          : (3,   None,       True,   False),
        'z_skewer'      : (4,   'Z',        False,  True)
    }

    @cached_property
    def id(self):
        return self.get_data((*self.data_dictionary['id']))
  

class PiccaStyleFiles(FilesSkewerBase):
    '''Class to handle PiccaStyle files output by LyaCoLoRe

    Attributes (extending FilesSkewerBase):
        lr_max (float): Maximum wavelength for masking.
        name (str): Name of the PiccaStyle file (e.g. picca-flux-noRSD-notnorm)
    '''

    data_dictionary = {
        'RA'            : (3,   'RA',       False,  False),
        'DEC'           : (3,   'DEC',      False,  False),
        'z'             : (3,   'Z',        False,  False),
        'id'            : (3,   'THING_ID', False,  False),
        'wavelength'    : (2,   None,       False,  True),
        'values'        : (0,   None,       True,   False)
    }

    def __init__(self,file_paths,lr_max, name, parent_sim, downsampling):
        '''Inits the PiccaStyleFiles object

        Args (extending FilesBase init args):
            lr_max (float): Maximum wavelength for masking.
            name (str): Name of the PiccaStyle file (e.g. picca-flux-noRSD-notnorm)
        '''
        FilesBase.__init__(self,file_paths, parent_sim, downsampling)
        self.__name__       = name 
        self.lr_max         = lr_max

    @cached_property
    def wavelength(self):
        return 10**self.get_data(*self.data_dictionary['wavelength'])

    @cached_property
    def z_skewer(self):
        return  (self.wavelength - utils.lya_rest)/utils.lya_rest
    
    @cached_property
    def mask(self):
        return  utils.make_IVAR_rows(self.lr_max,self.z,np.log10(self.wavelength))

    @cached_property
    def values(self):
        return self.get_data(*self.data_dictionary['values'])

class Spectra(FilesSkewerBase):
    '''Class to handle QuickQuasars output spectra files. 

    Attributes (extending FilesSkewerBase):
        arm (str): Arm where we are extracting the spectrum 
        spectra_files (list of str): Overriding file_paths from FilesBase
        truth_files (list of str): List of paths to the different truth_files
        zbest_files (list of str): List of paths to the different zbest_files
        lr_max (float): Set maximum wavelength for masking purposes 
        redhisft_to_use (str): We need a redshift for the Quasar to compute quantities (truth or best).

        truth (Files object): Object with the informatoin from the truth files.
        zbest (Files object): Object with the information from the zbest files. 

    '''

    data_dictionary = {
        'RA'            : (1,   'TARGET_RA',    False,  False),
        'DEC'           : (1,   'TARGET_DEC',   False,  False),
        'id'            : (1,   'TARGETID',     False,  False),
        'values'        : (0,   None,           True,   False)
    }

    def __init__(self, arm, spectra_files, truth_files, zbest_files, lr_max=1200., redshift_to_use='best', parent_sim=None, downsampling=1):
        '''Inits the Spectra object.

        Args (extending FilesBase init args):
            arm (str): Arm where we are extracting the spectrum 
            spectra_files (list of str): Overriding file_paths from FilesBase
            truth_files (list of str): List of paths to the different truth_files
            zbest_files (list of str): List of paths to the different zbest_files
            lr_max (float,optional): Set maximum wavelength for masking purposes 
            redhisft_to_use (str,optional): We need a redshift for the Quasar to compute quantities (truth or best).

        '''
        FilesBase.__init__(self,spectra_files, parent_sim, downsampling)
        self.truth      = TruthFiles(truth_files,  parent_sim)
        self.zbest      = BestFiles(zbest_files,   parent_sim)

        assert arm.lower() in ['r','b','z']
        assert redshift_to_use.lower() in ['truth','best']
        self.arm             = arm.lower()
        self.lr_max          = lr_max
        self.redshift_to_use = redshift_to_use.lower()

        self.data_dictionary['wavelength'] = ('{}_WAVELENGTH'.format(self.arm.upper()), None, False, True)
        self.data_dictionary['flux'] = ('{}_FLUX'.format(self.arm.upper()), None, True, False)

    @cached_property
    def z(self):
        if self.redshift_to_use == 'best':
            return self.zbest.z
        else:
            return self.truth.z

    @cached_property
    def z_skewer(self):
        return  (self.wavelength - utils.lya_rest)/utils.lya_rest
    
    @cached_property
    def mask(self):
        return  utils.make_IVAR_rows( self.lr_max, self.z, np.log10(self.wavelength) )

    @cached_property
    def flux(self):
        return self.get_data(*self.data_dictionary['flux'])