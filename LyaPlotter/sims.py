'''
    Module build to structure the information from different simulations. 
    
    The objective of this module is to allow an easier way of handling different simulations (of the same or different type).
'''

from lyacolore import utils
from LyaPlotter.file_types import *
import glob
import os
import re
import logging
from pathlib import Path
log = logging.getLogger(__name__)

def check_is_list(x): #pragma: no cover
    if isinstance(x,list) or isinstance(x,tuple):
        return True
    else:
        type_ = type(x)
        raise ValueError('A list/tuple was expected. Type of {} was {}'.format(x,type_),x)

class QuickQuasarsSim():
    '''A class used to handle QuickQuasars' output

    Attributes:
        id_ (int): Value to identify the specific simulation when dealing with several
        sim_path (str): Path to the QuickQuasars output
        nside (int): nside used in QuickQuasars
        __name__ (str): Name to identify the specific simulation when dealing with several
        file_type (str): Whether the simulation was a "gaussian" or a "2LPT" simulation
        compression (bool): Whether the ouptut files are compressed or not
        sim_class (str): Description of the class of simulation (set to "QuickQuasars"). 
    
    '''

    sim_class = 'QuickQuasars'

    def __init__(self,id_, sim_path, nside=16, name=None, file_type = 'gaussian', compression =True):
        """Inits QuickQuasarsSim without loading any information from the output.

        Args:
            id_ (int): Value to identify the specific simulation when dealing with several
            sim_path (str): Path to the QuickQuasars output
            nside (int, optional): nside used in QuickQuasars
            name (str, optional): Name to identify the specific simulation when dealing with several
            file_type (str, optional): Whether the simulation was a "gaussian" or a "2LPT" simulation
            compression (bool, optional): Whether the ouptut files are compressed or not.
        """
        self.id_        = id_
        self.__name__   = name if name else sim_path[sim_path.rfind('/')+1:]
        self.sim_path   = sim_path
        self.file_type  = file_type
        self.nside      = nside
        self.compression= compression

    def __str__(self): #pragma: no cover
        return "{} sim. Id: {}\tName: {}\tPath: {}".format(self.sim_class,self.id_,self.__name__,self.sim_path)

    @classmethod
    def search(cls,path):
        """Searches for QuickQuasars outputs.

        The search is performed by globbing a file with the format "/**/spectra-*-*.fits" which is found in the QuickQuasars output. This method has poor performance and shouldn't be used for large directories.

        Args:
            path (str): Path where to search for QuickQuasars outputs.

        Returns:
            A list with paths (str) to the different QuickQuasars outputs found. 
        """
        sim_paths = set()
        for path in glob.glob(path + '/**/spectra-*-*.fits',recursive=True):
            sim_paths.add( str(Path(path).parents[3]) )
        return sorted(list(sim_paths), key=lambda folder: os.stat(folder).st_ctime)

    def get_spectra(self, arm, pixels=[0], lr_max=1200., redshift_to_use= 'best', downsampling=1):
        ''' Get spectra for a single arm.

        It will load the different information from the fits files as they are called. 

        Args:
            arm (str): Arm for which we want to get spectra (r,b or z)
            pixels (list, optional): Pixels we want to use
            lr_max (float, optional): Set maximum wavelength for masking purposes 
            redhisft_to_use (str, optional): We need a redshift for the Quasar to compute quantities (truth or best).

        Returns:
            Spectra object for the given pixels.
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

        return Spectra(arm, spectra_files, truth_files, zbest_files, lr_max, redshift_to_use, self, downsampling)      

class CoLoReSim():
    ''' 
    A class used to handle CoLoRe output

    Attributes:
        id_ (int): Value to identify the specific simulation when dealing with several
        sim_path (str): Path to the CoLoRe output
        __name__ (str): Name to identify the specific simulation when dealing with several
        file_type (str): Whether the simulation was a "gaussian" or a "2LPT" simulation
        compression (bool): Whether the ouptut files are compressed or not
        sim_class (str): Description of the class of simulation (set to "CoLoRe"). 
    '''
    sim_class = 'CoLoRe'

    def __init__(self,id_,sim_path,name=None, file_type='gaussian'):
        '''Inits CoLoReSim without loading any information from the output.

        Args:
            id_ (int): Value to identify the specific simulation when dealing with several.
            sim_path (str): Path to the CoLoRe output.
            name (str, optional): Name to identify the specific simulation when dealing with several.
            file_type (str, optional): Whether the simulation is a "gaussian" or a "2LPT" simulation.
            
        '''
        self.id_ = id_

        self.sim_path = Path(sim_path)
        
        if name: #pragma: no cover
            self.__name__ = name
        else:
            self.__name__ = self.sim_path.name

        try:
            self.param_file = sorted( self.sim_path.glob('*cfg'), key=lambda x: x.stat().st_mtime )
        except:
            self.param_file = None
        self.file_type = file_type    

    @staticmethod
    def search(path):
        """Searches for CoLoRe outputs.

        The search is performed by globbing a file with the format "/**/out_srcs_s1_0.fits" which is found in the CoLoRe output. This method has poor performance and shouldn't be used for large directories.

        Args:
            path (str): Path where to search for CoLoRe outputs.

        Returns:
            A list with paths (str) to the different CoLoRe outputs found. 
        """
        sim_paths = []
        for path in glob.glob(path+'/**/out_srcs_s1_0.fits',recursive=True):
            sim_paths.append( path[:path.rfind('/')] )
        return sorted(sim_paths, key=lambda folder: os.stat(folder).st_ctime)

    def __str__(self): #pragma: no cover
        return "{} sim. Id: {}\tName: {}\tPath: {}".format(self.sim_class,self.id_,self.__name__,self.sim_path)

    def get_Sources(self, ifiles=None, lr_max=1200.,source=1,downsampling=1):
        ''' Get sources from a CoLoRe simulation.

        It will load the different information from the fits files as they are called.

        Args:
            ifiles (list of int, optional): Array with the srcs files that we want to use (default: all the files will be selected). 
            lr_max (float, optional): Maximum wavelength for masking.
            source (int, optional): Sources to analyse from the CoLoRe output.
            downsampling (float): Downsampling to apply to the data.

        Returns:
            A CoLoReFiles object.
        '''
        if ifiles is not None: 
            check_is_list(ifiles)
        else:
            files = self.sim_path.glob(f'out_srcs_s{ source }*')
            num = max( [int(re.search("(\d+).fits$",file).group(1)) for file in files] )
            ifiles = range(num+1)

        files            = [self.sim_path / 'out_srcs_s{}_{}.fits'.format(source,ifile) for ifile in ifiles]

        return CoLoReFiles(files, lr_max, self, downsampling=downsampling)
        # self.CoLoReFiles = CoLoReFiles(files,lr_max,self.file_type)
        # return self.CoLoReFiles

class LyaCoLoReSim():
    '''
    A class used to handle LyaCoLoRe output

    Attributes:
        id_ (int): Value to identify the specific simulation when dealing with several
        sim_path (str): Path to the LyaCoLoRe output
        nside (int): nside used in LyaCoLoRe
        __name__ (str): Name to identify the specific simulation when dealing with several
        file_type (str): Whether the simulation was a "gaussian" or a "2LPT" simulation
        compression (bool): Whether the ouptut files are compressed or not
        sim_class (str): Description of the class of simulation (set to "CoLoRe"). 

    '''
    sim_class = 'LyaCoLoRe'
    def __init__(self, id_, sim_path, nside=16, name=None, file_type='gaussian', compression=True):
        '''Inits LyaCoLoreSim without loading any information from the output.

        Args:
            id_ (int): Value to identify the specific simulation when dealing with several
            sim_path (str): Path to the LyaCoLoRe output
            nside (int, optional): nside used in LyaCoLoRe
            name (str, optional): Name to identify the specific simulation when dealing with several
            file_type (str, optional): Whether the simulation was a "gaussian" or a "2LPT" simulation
            compression (bool, optional): Whether the ouptut files are compressed or not
        '''
        self.id_        = id_
        self.__name__   = name if name else sim_path[sim_path.rfind('/')+1:]
        self.sim_path   = sim_path
        self.file_type  = file_type
        self.nside      = nside
        self.compression= compression
        self.picca_files= {}

    def __str__(self): #pragma: no cover
        return "{} sim. Id: {}\tName: {}\tPath: {}".format(self.sim_class,self.id_,self.__name__,self.sim_path)

    @classmethod
    def search(cls,path):
        """Searches for LyaCoLoRe outputs.

        The search is performed by globbing a file with the format "/**/master.fits" which is found in the CoLoRe output. This method has poor performance and shouldn't be used for large directories.

        Args:
            path (str): Path where to search for LyaCoLoRe outputs.

        Returns:
            A list with paths (str) to the different LyaCoLoRe outputs found. 
        """
        sim_paths = set()
        for path in glob.glob(path + '/**/master.fits',recursive=True):
            sim_paths.add( str(Path(path).parents[0]) )
        return sorted(list(sim_paths), key=lambda folder: os.stat(folder).st_ctime)

    def get_Transmission(self,  pixels=[0], lr_max=1200., downsampling=1):
        ''' Get transmission for the given pixels:

        Args:
            pixels (list of int, optional): Pixels we want to use
            lr_max (float, optional): Set maximum wavelength for mask
            downsampling (float): Downsampling to apply to the data.
        
        Returns:
            A TransmissionFiles object.
        '''
        check_is_list(pixels)
        files = []
        for pixel in pixels:
            dirname = utils.get_dir_name(self.sim_path, pixel)
            files.append( utils.get_file_name(dirname, 'transmission',     self.nside, pixel, self.compression) )

        return TransmissionFiles(files, lr_max, self, downsampling=downsampling)
        # self.Transmission = Transmission(files, lr_max, self.file_type)
        # return self.Transmission

    def get_GaussianCoLoRe(self,pixels=[0], lr_max=1200., downsampling=1):
        ''' Get colore skewers/velocity files for the given pixels:

        Args:
            pixels (list of int, optional): Pixels we want to use
            lr_max (float, optional): Set maximum wavelength for mask
            downsampling (float): Downsampling to apply to the data.

        Returns:
            A GaussianCoLoReFiles object
        '''
        check_is_list(pixels)
        files = []
        for pixel in pixels:
            dirname = utils.get_dir_name(self.sim_path, pixel)
            files.append(  utils.get_file_name(dirname, 'gaussian-colore', self.nside, pixel, self.compression)   )
        return GaussianCoLoReFiles(files, lr_max, self, downsampling=downsampling)
        # self.GaussianCoLoRe = GaussianCoLoRe(files, lr_max, self.file_type)
        # return self.GaussianCoLoRe

    def get_PiccaStyleFiles(self,file_name, pixels=[0], lr_max=1200, downsampling=1):
        ''' Get values stored in file file_name. It should be the format 'picca' files that we can get from LyaCoLoRe:

        Args:
            file_name (str): Name of the file (e.g. picca-flux-noRSD-notnorm)
            pixels (list of int, optional): Pixels we want to use
            lr_max (float, optional): Set maximum wavelength for mask
            downsampling (float): Downsampling to apply to the data.

        Returns:
            A PiccaStyleFiles object
        '''
        check_is_list(pixels)
        files = []
        for pixel in pixels:
            dirname = utils.get_dir_name(self.sim_path, pixel)
            files.append( utils.get_file_name(dirname, file_name,          self.nside, pixel, self.compression) )

        return PiccaStyleFiles(files, lr_max, file_name, self,downsampling=downsampling)
        # self.picca_files[file_name] = PiccaStyleFiles(files, lr_max,self.file_type,file_name)
        # return self.picca_files[file_name]