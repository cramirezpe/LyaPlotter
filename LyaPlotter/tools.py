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

def write_picca_drq_cat(Z, RA, DEC, out_file, THING_ID=None):

    from astropy.io import fits
    logger.info('Generating picca-like catalog from input data')
    logger.debug(f'Length of Z:\t{len(Z)}\nLength of RA:\t{len(RA)}\nLength of DEC:\t{len(DEC)}')
    logger.debug(f'Negative RA?\t{(RA<0).sum()}')
    THING_ID = np.arange(1, len(RA)+1) if THING_ID is None else THING_ID
    PLATE = THING_ID
    MJD = THING_ID
    FIBERID = THING_ID

    logger.info('Building fits columns')
    cols = []
    for value, label, dtype in zip((RA, DEC, Z, THING_ID, PLATE, MJD, FIBERID), ('RA', 'DEC', 'Z', 'THING_ID', 'MJD', 'FIBERID', 'PLATE'), ('D', 'D', 'D', 'K', 'K', 'K', 'K')):
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

def colore_to_drq_cat(in_sim, out_file, ifiles=None, source=1, downsampling=1, rsd=False, valid_pixels=None):
    ''' Method to extract a picca-like catalog from a CoLoRe box

        Args:
            in_path (str): Path to the CoLoRe box
            out_file (str): Path where to save the catalogue
            source (int, optional): Sources to analyse from the CoLoRe output.
            downsampling (float, optional): Downsampling to use when copying objects from one catalog to the other. (default: 1)
            rsd (bool): Apply rsd
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
    write_picca_drq_cat(Z, RA, DEC, out_file)

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
    THING_ID = h[1]['THING_ID'][:]

    pixels = healpy.ang2pix(nside, RA, DEC, lonlat=True)

    for i in range(12*nside**2):
        logger.info(f'Computing catalog for pixel {i}')
        file = Path(out_path) / f'pixel_{i}.fits'
        mask = pixels == i
        logger.info(f'\tLength of catalog: {mask.sum()}')
        write_picca_drq_cat(Z[mask], RA[mask], DEC[mask], out_file=file, THING_ID=THING_ID[mask])


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