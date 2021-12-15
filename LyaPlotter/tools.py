'''
    Module with different tools not directly related with plots. 
'''

import argparse
import fitsio
import scipy as sp
import healpy
import numpy as np
import os
from pathlib import Path
import logging
logger = logging.getLogger(__name__)

def write_picca_drq_cat(Z, RA, DEC, out_file, THING_ID=None, simplified=False):

    from astropy.io import fits
    logger.info('Generating picca-like catalog from input data')
    logger.debug(f'Length of Z:\t{len(Z)}\nLength of RA:\t{len(RA)}\nLength of DEC:\t{len(DEC)}')
    logger.debug(f'Negative RA?\t{(RA<0).sum()}')

    if not simplified:
        THING_ID = np.arange(1, len(RA)+1) if THING_ID is None else THING_ID
        PLATE = THING_ID
        MJD = THING_ID
        FIBERID = THING_ID

    logger.info('Building fits columns')
    cols = []
    
    if simplified:
        values = (RA,DEC,Z)
        labels = ('RA', 'DEC', 'Z')
        dtypes = ('D', 'D', 'D')
    else:
        values = (RA, DEC, Z, THING_ID, PLATE, MJD, FIBERID)
        labels = ('RA', 'DEC', 'Z', 'THING_ID', 'MJD', 'FIBERID', 'PLATE')
        dtypes = ('D', 'D', 'D', 'K', 'K', 'K', 'K')

    for value, label, dtype in zip(values, labels, dtypes):
        cols.append(fits.Column(name=label, format=dtype, array=value))

    logger.debug('Defining FITS headers')
    p_hdr = fits.Header()
    t_hdr = fits.Header()

    t_hdr['NSIDE'] = 16

    logger.debug('Defining hdulist')
    hdulist = fits.HDUList()

    hdulist.append(fits.PrimaryHDU(header=p_hdr))
    hdulist.append(fits.BinTableHDU.from_columns(cols, header=t_hdr))

    logger.info('Writting data to file')
    hdulist.writeto(out_file, overwrite=True)
    hdulist.close() 

def colore_to_drq_cat(in_sim, out_file, ifiles=None, source=1, downsampling=1, rsd=False, valid_pixels=None, simplified=False):
    ''' Method to extract a picca-like catalog from a CoLoRe box

        Args:
            in_path (str): Path to the CoLoRe box
            out_file (str): Path where to save the catalogue
            source (int, optional): Sources to analyse from the CoLoRe output.
            downsampling (float, optional): Downsampling to use when copying objects from one catalog to the other. (default: 1)
            rsd (bool): Apply rsd
            simplified (bool, optional): Use simplified catalogs (only RA DEC z)
    '''
    from LyaPlotter.sims import CoLoReSim
    from astropy.io import fits

    in_sim = Path(in_sim)
    sim = CoLoReSim(0, in_sim)
    data = sim.get_Sources(ifiles=ifiles, source=source, downsampling=downsampling)

    if rsd:
        Z = data.z + data.dz_rsd
    else:
        Z = data.z
    RA = data.RA
    DEC = data.DEC

    if valid_pixels is not None:
        print('valid_pixels', valid_pixels)
        mask = np.in1d(data.healpixels, np.asarray(valid_pixels))
        print('objects in valid pixels:', mask.sum())
        Z = Z[mask]
        RA = data.RA[mask]
        DEC = data.DEC[mask]

    print(f'Length of resulting catalog: \t{len(Z)}')
    write_picca_drq_cat(Z, RA, DEC, out_file, simplified=simplified)

def trim_catalog_into_pixels(in_cat, out_path, nside):
    ''' Method to trim a catalog into multiple catalogs (one for each healpix pixel)

    Args:
        in_cat (str or Path): Input catalog
        out_path (str or Path): Output path for subcatalogs
        nside (int): pixelization that will determine the number of subcatalogs (npix = 12*nside**2)
    '''

    Path(out_path).mkdir(parents=True, exist_ok=False)

    logger.info('Reading input catalog')
    h = fitsio.FITS(in_cat)

    logger.info(f'Length of cat: {h[1].get_nrows()}')
    
    Z = h[1]['Z'][:]
    RA = h[1]['RA'][:]
    DEC = h[1]['DEC'][:]


    pixels = healpy.ang2pix(nside, RA, DEC, lonlat=True)

    for i in range(12*nside**2):
        logger.info(f'Computing catalog for pixel {i}')
        file = Path(out_path) / f'pixel_{i}.fits'
        mask = pixels == i
        logger.info(f'\tLength of catalog: {mask.sum()}')
        write_picca_drq_cat(Z[mask], RA[mask], DEC[mask], out_file=file, THING_ID=None, simplified=True)

def generate_random_objects_from_nz_file(dndz_file, out_file, zmin=None, zmax=None, N_factor=1, pixel_footprint=None, nside=2, method='full_sky_dice'):
    ''' Generate random objects from a given dndz distribution.

    Args:
        dndz_file (Path): Input filename to read dndz.
        out_file (Path): Output path to save the random catalog.
        zmin (float, optional): Min. redshift. (Default: min redshift in file).
        zmax (float, optional): Max. redshift. (Default: max redshift in file).
        N_factor (float, optional): N_rand/N_data assuming N_data from dndz file, and dndz in units of deg^-2.  (Default: 1)
        pixel_footprint (array of int, optional): Footprint randoms should span. (Default: full sky)
        nside (int, optional): Nside for the previous footprint pixels. (Default: 2)
        method (str, optonal): Method to generate random positions in the sky. Options:
            - full_sky_dice (throw randoms at the whole sky and then removing the ones not matching). (Default)
            - pixel (make random pixel by pixel). 
    '''
    from scipy.interpolate import interp1d 
    from scipy.integrate import quad

    dndz_file = Path(dndz_file)

    logger.info(f'Generating random catalog from file: {str(dndz_file)}')
    in_z, in_nz = np.loadtxt(dndz_file, unpack=True)

    zmin = zmin if zmin != None else in_z[0]
    zmax = zmax if zmax != None else in_z[-1]

    # Cumulative distribution
    in_Nz = np.cumsum(in_nz)
    N = interp1d(in_z, in_Nz)
    N_inv = interp1d(in_Nz, in_z)
    
    pixarea = healpy.pixelfunc.nside2pixarea(nside, degrees=True)
    if pixel_footprint is None:
        pixels = 48
    else:
        pixels = len(pixel_footprint)
    
    area = pixarea*pixels
    interpolation = interp1d(in_z, in_nz)
    NRAND = int(quad(interpolation, zmin, zmax)[0]*area)

    logger.debug('Generating random redshifts')
    ran = np.random.random(NRAND)
    z = N_inv( ran*(N(zmax)-N(zmin)) + N(zmin) )

    RA, DEC = generate_random_positions(NRAND, pixels=pixel_footprint, nside=nside, method=method)

    write_picca_drq_cat(z, RA, DEC, out_file, simplified=True)

def generate_random_objects_from_data(z, RA, DEC, out_file, N_factor=1, pixel_footprint=None, nside=2, method='full_sky_dice'):
    ''' Generate random objects from a given catalogue
    
    Args:
        z (array of float): Redshift values from catalog.
        RA (array of float): Right Ascention values from catalog.
        DEC (array of float): Declination values from catalog.
        out_file (Path): Output path to save the random catalog.
        N_factor (float, optional): N_rand/N_data assuming N_data from dndz file, and dndz in units of deg^-2.  (Default: 1)
        pixel_footprint (array of int, optional): Footprint randoms should span. (Default: full sky)
        nside (int, optional): Nside for the previous footprint pixels. (Default: 2)
        method (str, optonal): Method to generate random positions in the sky. Options:
            - full_sky_dice (throw randoms at the whole sky and then removing the ones not matching). (Default)
            - pixel (make random pixel by pixel). 
    '''
    from scipy.interpolate import interp1d

    logger.info('Generating random catalog')
    logger.info(f'Sorting by redshift')
    
    z_sort = np.sort(z)
    NDATA = len(z_sort)
    NRAND = int(NDATA*N_factor)
    
    logger.info(f'Data length: {NDATA}')
    logger.info(f'Randoms length: {NRAND}')
    
    logger.info('Creating inverse PDF')
    p = np.linspace(0, 1, NDATA, endpoint=True) #endpoint set to true will cause a biased estimator... but I chose it anyway to avoid invalid values later on.

    z_gen = interp1d(p, z_sort)

    logger.info('Interpolating redshift')
    ran1 = np.random.random(int(NRAND))
    z = z_gen(ran1)
    
    RA, DEC = generate_random_positions(NRAND, pixels=pixel_footprint, nside=nside, method=method)
    
    write_picca_drq_cat(z, RA, DEC, out_file, simplified=True)

def generate_random_positions(NRAND, pixels=None, nside=2, method='full_sky_dice'):
    '''Generate random positions in the sky
    
    Args:
        NRAND (int): Number of positions to generate.
        pixels (array of int, optional): Pixels to generate randoms in. (Default: None -> all sky).
        nside (int, optional): Nside for the pixelization of pixels. (Default: 2).
        method (str, optional): Method to generate random positions. Options:
            - full_sky_dice: throw randoms at the whole sky and then remove the ones not matching any pixel. (Default).
            - pixel: Generate randoms around each pixel, discarding less number of randoms if the footprint is small enough.
    ''' 
    logger.info(f'Computing random positions in the sky')
    if pixels is None:
        logger.info('No pixel mask added. Computing pixels for the full sky')
        ran1 = np.random.random(int(NRAND))
        ran2 = np.random.random(int(NRAND))

        RA = np.degrees(np.pi*2*ran1)
        DEC = np.degrees(np.arcsin(2.*(ran2-0.5)))

        return RA, DEC
    else:
        if method == 'full_sky_dice':
            logger.info('Computing randoms for full sky and discarding the ones not in footprint')
            sky_fraction = len(pixels) / (12*nside**2)
            logger.debug(f'Sky fraction: {sky_fraction}')

            RA = []
            DEC = []

            while len(RA) < NRAND:
                logger.debug('Iterating through the random generator')
                randoms_left = NRAND - len(RA)
                randoms_to_generate = int(randoms_left/sky_fraction)
                logger.debug(f'\t{randoms_left} elements left to include\n\tGenerating {randoms_to_generate} random positions')
                ran1 = np.random.random(randoms_to_generate)
                ran2 = np.random.random(randoms_to_generate)

                RAs     = np.degrees(np.pi*2*ran1)
                DECs    = np.degrees(np.arcsin(2.*(ran2-0.5)))

                gen_pixels = healpy.pixelfunc.ang2pix(nside, RAs, DECs, lonlat=True)
                w = np.isin(gen_pixels, pixels)
                logger.debug(f'Accepted randoms: {w.sum()}')

                RA = np.append(RA, RAs[w])
                DEC = np.append(DEC, DECs[w])
            
            return RA, DEC
        
        elif method == 'pixel':
            logger.info('Computing random positions around valid pixels')

            # First thing is to poisson sampling the pixels I have
            # I'll need to match exactly the number of input pixels
            # for convenience, therefore I'd need to add extra obj
            # until the overall number of randoms is matched
            _lambda = NRAND / len(pixels)
            randoms_per_pixel = np.random.poisson(_lambda, len(pixels))
            extra_objs = randoms_per_pixel.sum() - NRAND
            if extra_objs > 0:
                for i in range(np.abs(extra_objs)):
                    randoms_per_pixel[np.random.randint(0, len(pixels))] -= 1
            elif extra_objs < 0:
                for i in range(np.abs(extra_objs)):
                    randoms_per_pixel[np.random.randint(0, len(pixels))] += 1

            pix_area = healpy.pixelfunc.nside2pixarea(nside)
            _index = 0
            RA = []
            DEC = []

            for pixel, N in zip(pixels, randoms_per_pixel):
                logger.info(f'Computing random positions for pixel {pixel}')
                pixel_center = healpy.pix2ang(nside=nside, ipix=pixel)
                corners = healpy.rotator.vec2dir(
                    healpy.boundaries(nside, [pixel], step=1000)[0]
                )
                th_range = (min(corners[0]), max(corners[0]))
                ph_range = (min(corners[1]), max(corners[1]))

                range_size = np.abs((ph_range[1]-ph_range[0]) * (np.cos(th_range[0])-np.cos(th_range[1])))

                TH = []
                PH = []
                valid_fraction = pix_area / range_size
                while len(TH) < N:
                    randoms_left = N - len(PH)
                    logger.debug(f'Randoms left: {randoms_left}')
                    ran1 = np.random.random(int(randoms_left/valid_fraction))
                    ran2 = np.random.random(int(randoms_left/valid_fraction))

                    PHs = ran1*(ph_range[1]-ph_range[0]) + ph_range[0]
                    THs = np.arccos( np.cos(th_range[0]) - ran2*( np.cos(th_range[0])-np.cos(th_range[1]) ) )
                    
                    new_pixels = healpy.pixelfunc.ang2pix(nside, THs, PHs)
                    mask = new_pixels == pixel
                    TH = np.append(TH, THs[mask])
                    PH = np.append(PH, PHs[mask])

                RA = np.append(RA, np.degrees(PH))
                DEC = np.append(DEC, 90 - np.degrees(TH))
            
            return RA, DEC

        else:
            raise ValueError('Select a valid method.', method)


def master_to_qso_cat(in_file, out_file, min_zcat=1.7, nside=16, downsampling=1, downsampling_seed=0, randoms_input=False, rsd=True):
    ''' 
        Method to convert master.fits files from LyaCoLoRe sims into qso catalogs zcat.fits.

        Args:
            in_file (str): Path of master.fits file used as input.
            out_file (str): Path where to save the z_cat.fits output.
            min_zcat (float, optional): Min z of the objects that will be copied from one catalog to the other. (default: 1.7)
            nside (int, optional): Nside used in the LyaCoLoRe simulation. (default: 16)
            downsampling (float, optional): Downsampling to use when copying objects from one catalog to the other. (default: 1)
            downsampling_seed(int,optional): (default: 0)
            randoms_input (bool,optional): Whether the input catalog is a catalog of randoms or not. Basically to select 'Z' as redshift instead of 'Z_QSO_RSD'. (default: False)

        Returns:
            None
    '''
    if randoms_input:
        z = 'Z'
    elif rsd:
        z = 'Z_QSO_RSD'
    else:
        z = 'Z_QSO_NO_RSD'

    ### Make random generator
    state = sp.random.RandomState(downsampling_seed)

    ### Data
    h = fitsio.FITS(in_file)
    m_data = sp.sort(h[1].read(),order=['MOCKID',z])
    data = {}
    for k in ['RA','DEC']:
        data[k] = m_data[k][:]
    for k in ['THING_ID','PLATE','MJD','FIBERID']:
        data[k] = m_data['MOCKID'][:]
    data['Z'] = m_data[z][:]
    print(data['Z'].min())
    w = data['Z']>min_zcat
    for k in data.keys():
        data[k] = data[k][w]
    h.close()
    phi = data['RA']*sp.pi/180.
    th = sp.pi/2.-data['DEC']*sp.pi/180.
    pix = healpy.ang2pix(nside,th,phi)
    data['PIX'] = pix
    print('INFO: {} QSO in mocks data'.format(data['RA'].size))

    ### Get reduced data numbers
    original_nbData = data['RA'].shape[0]
    nbData = round(original_nbData * downsampling)

    ### Save data
    assert nbData<=data['RA'].size
    w = state.choice(sp.arange(data['RA'].size), size=nbData, replace=False)
    w_thid = data['THING_ID'][w]
    print(w_thid.shape)
    print('INFO: downsampling to {} QSOs in catalog'.format(nbData))
    out = fitsio.FITS(out_file,'rw',clobber=True)
    cols = [ v[w] for k,v in data.items() if k not in ['PIX'] ]
    names = [ k for k in data.keys() if k not in ['PIX'] ]
    out.write(cols,names=names)
    out.close()

    return

def from_string_to_options(input):
    '''
        Function to convert file or string with argparse options into a dict of these options.

        Args: 
            input (str): Input file or string which contains the argparse options (e.g: --option1 param_a --option2 param_b --store_true_option)

        Returns:
            (dict): Dictionary containing the different options. 
    '''

    if os.path.isfile(input):
        with open(input, 'r') as file:
            string = file.read().strip()
    else:
        string = input.strip()

    options = dict()
    for x in string.split('--')[1:]:
        x = x.strip().split(' ',1)
        try:
            options[x[0]] = x[1]
        except IndexError:
            options[x[0]] = ""
    return options