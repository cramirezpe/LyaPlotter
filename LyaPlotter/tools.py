'''
    Module with different tools not directly related with plots. 
'''

import fitsio
import scipy as sp
import healpy
import os

def master_to_qso_cat(in_file, out_file, min_zcat=1.7, nside=16, downsampling=1, downsampling_seed=0, randoms_input=False):
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
    z = 'Z' if randoms_input else 'Z_QSO_RSD'

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