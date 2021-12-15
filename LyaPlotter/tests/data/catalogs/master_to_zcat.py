''' This is a copypaste from James code to create a test version'''

import fitsio
import scipy as sp
import healpy


class empty():
    pass

args = empty()

args.downsampling_seed = 2
args.in_dir = '/global/homes/c/cramirez/Codes/LyaPlotter/LyaPlotter/tests/data/catalogs'
args.min_cat_z = 0
args.nside = 16
args.out_dir = '/global/homes/c/cramirez/Codes/LyaPlotter/LyaPlotter/tests/data/catalogs'
args.downsampling = 0.5

def create_qso_cat(args):

    ### Make random generator
    state = sp.random.RandomState(args.downsampling_seed)

    ### Data
    h = fitsio.FITS(args.in_dir+'/master.fits')
    m_data = sp.sort(h[1].read(),order=['MOCKID','Z_QSO_RSD'])
    data = {}
    for k in ['RA','DEC']:
        data[k] = m_data[k][:]
    for k in ['THING_ID','PLATE','MJD','FIBERID']:
        data[k] = m_data['MOCKID'][:]
    data['Z'] = m_data['Z_QSO_RSD'][:]
    print(data['Z'].min())
    w = data['Z']>args.min_cat_z
    for k in data.keys():
        data[k] = data[k][w]
    h.close()
    phi = data['RA']*sp.pi/180.
    th = sp.pi/2.-data['DEC']*sp.pi/180.
    pix = healpy.ang2pix(args.nside,th,phi)
    data['PIX'] = pix
    print('INFO: {} QSO in mocks data'.format(data['RA'].size))

    ### Get reduced data numbers
    original_nbData = data['RA'].shape[0]
    nbData = round(original_nbData * args.downsampling)

    ### Save data
    assert nbData<=data['RA'].size
    w = state.choice(sp.arange(data['RA'].size), size=nbData, replace=False)
    w_thid = data['THING_ID'][w]
    print(w_thid.shape)
    print('INFO: downsampling to {} QSOs in catalog'.format(nbData))
    out = fitsio.FITS(args.out_dir+'/zcat_{}.fits'.format(args.downsampling),'rw',clobber=True)
    cols = [ v[w] for k,v in data.items() if k not in ['PIX'] ]
    names = [ k for k in data.keys() if k not in ['PIX'] ]
    out.write(cols,names=names)
    out.close()

    return

create_qso_cat(args)