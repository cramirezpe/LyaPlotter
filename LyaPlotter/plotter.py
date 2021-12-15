'''
    Module build to help making plots
'''

import matplotlib.pyplot as plt
import numpy as np
from LyaPlotter.computations import Computations
import healpy as hp

class Plotter: #pragma: no cover
    @staticmethod
    def plot_footprint(RA,DEC,bins,ax=None,**kwargs):
        if not ax: fig, ax = plt.subplots()

        my_cmap = plt.cm.jet
        my_cmap.set_under('w',1)

        plt.hist2d(RA,DEC,bins=bins, vmin=1, cmap=my_cmap)
        cb = plt.colorbar()
        cb.set_label('Number of entries')
        return
    
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
    def plot_dist(values,bins,ax=None, **kwargs):
        if not ax: fig, ax = plt.subplots()
        ax.hist(values,bins, **kwargs)
        return

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

    def plot_mollview_positions(RA, DEC, nside=64, ax=None, **kwargs):
        '''
            Plot mollview of objects positions

            Args:
                RA: (in degrees)
                DEC: (in degrees)
                nside:
                ax: Set to None to create new axis.
                kwargs: kwargs to be sent to mollview
        '''
        if ax is not None:
            plt.axes(ax)
            hold = True
        else:
            hold = False
        
        npix = hp.nside2npix(nside)
        pix = hp.ang2pix(nside, np.radians(90-DEC), np.radians(RA))
        n = np.bincount(pix, minlength=npix)
        
        hp.mollview(n, hold = hold)