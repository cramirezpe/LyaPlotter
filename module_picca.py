import scipy as sp
import fitsio
import matplotlib.pyplot as plt
import picca.wedgize
import glob 
import os
from pathlib import Path

class PlotterSkeleton():
    def __init__(self,path,name=None):
        self.__name__   = name if name else path
        self.path       = path
        self.hdulist    = fitsio.FITS(path)

    @classmethod
    def search(cls, path, file_name):
        sim_paths = set()
        for path in glob.glob(path + '/**/'+file_name, recursive=True):
            sim_paths.add( path )
        return sorted(list(sim_paths), key=lambda folder: os.stat(folder).st_ctime)

class DeltasPlotter(PlotterSkeleton):
    ''' 
    Class to store the usual plots from the delta_attributes.fits.gz output. It should be initialized through:
        object = DeltasPlotter(path, name)
    If the name is not provided it would default to the folder's name. 

    To define the data needed for plotting it should be run:
        object.get_data()

    After that, the following plot methods are available (kwargs sent to the main plot):
        plot_stack(self, ax=None, **kwargs)
        plot_eta(self, ax=None, **kargs)
        plot_var_lss(self, ax=None, **kwargs)
        plot_fudge(self, ax=None, **kwargs)
        plot_mean_cont(self, ax=None, **kwargs)
    '''
    @classmethod
    def search(cls,path):
        ''' Classmethod to search delta_attributes outputs under the given path. They are returned sorted by date.''' 
        return PlotterSkeleton.search(path,'delta_attributes.fits.gz')

    def __str__():
        return 0
    
    def get_data(self):
        self.stack_loglam   = self.hdulist[1]['LOGLAM'][:]
        self.stack          = self.hdulist[1]['STACK'][:]
        self.stack_weight   = self.hdulist[1]['WEIGHT'][:]
        self.weight_loglam  = self.hdulist[2]['LOGLAM'][:]
        self.weight_eta     = self.hdulist[2]['ETA'][:]
        self.weight_nb_pixels=self.hdulist[2]['NB_PIXELS'][:]
        self.weight_var_lss = self.hdulist[2]['VAR_LSS'][:]
        self.weight_fudge   = self.hdulist[2]['FUDGE'][:]
        self.cont_mean      = self.hdulist[3]['MEAN_CONT'][:]
        self.cont_loglam_rest=self.hdulist[3]['LOGLAM_REST'][:]
        self.cont_weight    = self.hdulist[3]['WEIGHT'][:]
    
        self.cut_stack  = (self.stack!=0.)           & (self.stack_weight >0.)
        self.cut_eta    = (self.weight_nb_pixels>0.) & (self.weight_eta != 1.)
        self.cut_var_lss= (self.weight_nb_pixels>0.) & (self.weight_var_lss!=0.1)
        self.cut_fudge  = (self.weight_nb_pixels>0.) & (self.weight_fudge!=1.e-7)
        self.cut_cont   = (self.cont_mean!=0.)       & (self.cont_weight>0.)
    
    def plot_stack(self, ax=None, **kwargs):
        if not ax: fig, ax = plt.subplots()
        loglam = self.stack_loglam[self.cut_stack]
        stack  = self.stack[self.cut_stack]
        ax.plot(10.**loglam,stack, **kwargs)
        ax.grid()
        ax.set_xlabel(r'$\lambda_{\mathrm{Obs.}} \, [\AA]$')
        ax.set_ylabel(r'$\mathrm{\overline{Flux}}$')
        return

    def plot_eta(self, ax=None, **kwargs):
        if not ax: fig, ax = plt.subplots()
        loglam = self.weight_loglam[self.cut_eta]
        eta    = self.weight_eta[self.cut_eta]
        ax.errorbar(10.**loglam, eta, **kwargs)
        ax.grid()
        ax.set_xlabel(r'$\lambda_{\mathrm{Obs.}} \, [\AA]$')
        ax.set_ylabel(r'$\eta$')
        return
    
    def plot_var_lss(self, ax=None, **kwargs):
        if not ax: fig, ax = plt.subplots()
        loglam = self.weight_loglam[self.cut_var_lss]
        var_lss= self.weight_var_lss[self.cut_var_lss]
        ax.errorbar(10.**loglam, var_lss, **kwargs)
        ax.grid()
        ax.set_xlabel(r'$\lambda_{\mathrm{Obs.}} \, [\AA]$')
        ax.set_ylabel(r'$\sigma^{2}_{\mathrm{LSS}}$')
        return

    def plot_fudge(self, ax=None, **kwargs):
        if not ax: fig, ax = plt.subplots()
        loglam = self.weight_loglam[self.cut_fudge]
        fudge  = self.weight_fudge[self.cut_fudge]
        ax.errorbar(10.**loglam, fudge, **kwargs)
        ax.grid()
        ax.set_xlabel(r'$\lambda_{\mathrm{Obs.}} \, [\AA]$')
        ax.set_ylabel(r'$\mathrm{Fudge}$')
        return

    def plot_mean_cont(self, ax=None, **kwargs):
        if not ax: fig, ax = plt.subplots()
        loglam  = self.cont_loglam_rest[self.cut_cont]
        mean_cont = self.cont_mean[self.cut_cont]
        ax.plot(10.**self.cont_loglam_rest, mean_cont,**kwargs)
        ax.grid()
        ax.set_xlabel(r'$\lambda_{\mathrm{R.F.}} \, [\AA]$')
        ax.set_ylabel(r'$\mathrm{\overline{Flux}}$')
        return 

class cfPlotter(PlotterSkeleton):
    ''' 
    Class to store the usual plots from the e_cf.fits.gz output. It should be initialized through:
        object = cfPlotter(path, name)
    If the name is not provided it would default to the folder's name. 

    After that, the following plot methods are available (kwargs sent to the main plot):
        plot(self, rfactor=0, mumin=0, mumax=0.5, ax=None, title=False, **kwargs)
    '''
    @classmethod
    def search(cls,path):
        ''' Classmethod to search outputs under the given path. They are returned sorted by date.''' 
        return PlotterSkeleton.search(path,'e_cf.fits.gz')

    def plot(self, rfactor=0, mumin=0., mumax=1, ax=None, title=False, wedge_args={}, **kwargs):
        da     = self.hdulist[1]['DA'][:]
        co     = self.hdulist[1]['CO'][:]

        if not ax: fig, ax = plt.subplots()
        w = picca.wedgize.wedge(mumin=mumin,mumax=mumax, **wedge_args) # To make plots of other wedges modify accordingly
        w.wedge(da,co)
        self.data_wedge = w.wedge(da,co)
        coef = self.data_wedge[0]**rfactor

        ax.errorbar(self.data_wedge[0], coef*self.data_wedge[1], yerr=coef*sp.sqrt(sp.diag(self.data_wedge[2])), **kwargs)        
        ax.grid()
        ax.set_xlabel(r"$r \, [\mathrm{Mpc \, h^{-1}}]$")
        ax.set_ylabel(r"$r^{0} \xi(r) \, [\mathrm{{Mpc \, h^{{-1}}  }}]$".format(rfactor))
        if title: ax.set_title('{0} < $\mu$ < {1}'.format(mumin, mumax))
        return

class cf1dPlotter(PlotterSkeleton):
    ''' 
    Class to store the usual plots from the cf.fits.gz output. It should be initialized through:
        object = cf1dPlotter(path, name)
    If the name is not provided it would default to the folder's name. 

    After that, the following plot methods are available (kwargs sent to the main plot):
        normalized_correlation(self, ax=None, **kwargs)
        evolution_variance(self, ax=None, **kwargs)
    '''
    def __init__(self,path,name=None):
        self.__name__ = name
        self.path = path
        self.data = fitsio.FITS(path)
        self.head  = self.data[1].read_header()
        self.dll   = self.head['DLL']
        
    @classmethod
    def search(cls,path):
        ''' Classmethod to search outputs under the given path. They are returned sorted by date.''' 
        return PlotterSkeleton.search(path,'cf.fits.gz')

    def evolution_variance(self, ax=None, **kwargs):
        if not ax: fig, ax = plt.subplots()
        llmin = self.head['LLMIN']
        llmax = self.head['LLMAX']

        n1d   = int((llmax-llmin)/self.dll+1)

        x    = sp.arange(0.,n1d)*self.dll+llmin
        y    = self.data[1]['v1d'][:]
        keep = (self.data[1]['wv1d'][:]>0.)
        x    = x[keep]
        y    = y[keep]
        ax.plot(10.**x,y,**kwargs)

        ax.set_xlabel(r'$\lambda_{\mathrm{Obs.}} \, [\AA] $')
        ax.set_ylabel(r'$\xi^{ff,1D}$')
        ax.grid()
        return
            
    def normalized_correlation(self, ax=None, **kwargs):
        if not ax: fig, ax = plt.subplots()
        y    = self.data[1]['c1d'][:]
        binsize = self.dll
        bins = sp.arange(y.size)
        x = sp.power(10,bins*binsize)
        keep = (self.data[1]['nc1d'][:]>0.)
        x = x[keep]
        y = y[keep]
        ax.plot(x,y,**kwargs)
        ax.set_xlabel(r'$\lambda_{1}/\lambda_{2}$')
        ax.set_ylabel(r'$\xi^{ff,1D}(1,2) / \sqrt(\xi^{ff,1D}(1)\xi^{ff,1D}(2))}$')
        ax.grid()
        return

class xcfPlotter(PlotterSkeleton):
    ''' 
    Class to store the usual plots from the e_xcf.fits.gz output. It should be initialized through:
        object = cf1dPlotter(path, name)
    If the name is not provided it would default to the folder's name. 

    After that, the following plot methods are available (kwargs sent to the main plot):
        plot(self, rfactor=0, mumin=0., mumax=1., ax=None, title=False, **kwargs)
    '''
    @classmethod
    def search(cls,path):
        ''' Classmethod to search outputs under the given path. They are returned sorted by date.''' 
        return PlotterSkeleton.search(path,'e_xcf.fits.gz')
      
    def plot(self, rfactor=0, mumin=0., mumax=1., ax=None, title=False, **kwargs):
        da = self.hdulist[1]['DA'][:]
        co = self.hdulist[1]['CO'][:]
        w = picca.wedgize.wedge(rpmin=-200.,rpmax=200.,nrp=100,rtmin=0.,rtmax=200.,nrt=50,mumin=mumin,mumax=mumax)
        w.wedge(da,co)
        data_wedge = w.wedge(da,co)
        self.data_wedge = data_wedge
        coef = data_wedge[0]**rfactor
        ax.errorbar(data_wedge[0],coef*data_wedge[1],yerr=coef*sp.sqrt(sp.diag(data_wedge[2])),**kwargs)
        ax.grid()
        ax.set_xlabel(r"$r \, [\mathrm{Mpc \, h^{-1}}]$")
        ax.set_ylabel(r"$r^{0} \xi(r) \, [\mathrm{{Mpc \, h^{{-1}}  }}]$".format(rfactor))
        if title: ax.set_title('{0} < $\mu$ < {1}'.format(mumin, mumax))
        return

class Searcher:
    def __init__(self,path,type_):
        types = ['deltas','cf1d','cf','xcf']
        if type_ not in types:
            raise ValueError('Type {} not recognized. Type should be in: {}'.format(type_,types))
        self.path = path
        self.type_ = type_
        
    def search(self,file_name=None):
        sims = {}
        
        if self.type_ == 'deltas':
            if not file_name: file_name = 'delta_attributes.fits.gz'
            Plotter = DeltasPlotter
        elif self.type_ == 'cf1d':
            if not file_name: file_name = 'Correlations/cf.fits.gz'
            Plotter = cf1dPlotter
        elif self.type_ == 'cf':
            if not file_name: file_name = 'e_cf.fits.gz'
            Plotter = cfPlotter
        elif self.type_ == 'xcf':
            if not file_name: file_name = 'e_xcf.fits.gz'
            Plotter = xcfPlotter
        
        for i,fits_file in enumerate(glob.glob(self.path+'/**/'+file_name, recursive=True)):
            name = fits_file[len(self.path):]
            name = name[:name.find('/')]
            print("Id: {}\tName: {}\tPath: {}".format(i,name,fits_file))
            sims[i] = Plotter(fits_file,name=name)
        self.sims = sims 