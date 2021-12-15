import argparse
from pathlib import Path
from LyaPlotter.tools import master_to_qso_cat, colore_to_drq_cat, trim_catalog_into_pixels, generate_random_objects_from_data, generate_random_objects_from_nz_file
import logging
import sys

def master_to_qso_cat_script(): #pragma: no cover
    parser = argparse.ArgumentParser(description='Convert master.fits file from LyaCoLoRe sim into zcat.fits file that can be read by picca')

    parser.add_argument('--master', required=True, help='Path to master.fits file from LyaCoLoRe sim')
    parser.add_argument('--out-file', required=True, help='Path to output zcat file')
    parser.add_argument('--z-min', default=1.7, type=float)
    parser.add_argument('--nside', default=16, type=int, help='nside used in the LyaCoLoRe box')
    parser.add_argument('--downsampling', type=float, default=1)
    parser.add_argument('--downsampling-seed', type=int, default=0)
    parser.add_argument('--rsd', action='store_true', help='Include RSD distortions in resulting catalog.')
    parser.add_argument('--randoms-input', action='store_true', help='Set this option if the input catalog is a randoms (Will select "Z" as redshift instead of "Z_QSO_RSD"')
    parser.add_argument('--log-level', default='WARNING', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'NOTSET'])

    args = parser.parse_args()
    level = logging.getLevelName(args.log_level)
    logging.basicConfig(stream=sys.stdout, level=level)

    master_to_qso_cat(in_file=args.master, out_file=args.out_file, min_zcat=args.z_min, nside=args.nside, downsampling=args.downsampling, downsampling_seed=args.downsampling_seed, randoms_input=args.randoms_input, rsd=args.rsd)

def colore_to_drq_script(): #pragma: no cover
    parser = argparse.ArgumentParser(description='Method to extract a picca-like catalog from a CoLoRe box')
    
    parser.add_argument('--in-sim',    required=True, type=Path, help='Path to the CoLoRe box.')
    parser.add_argument('--out-file',   required=True, type=Path, help='Path where to save the catalogue.')
    parser.add_argument('--source',     default=1, type=int, help='Source to analyse from the CoLoRe output.')
    parser.add_argument('--ifiles', default=None, nargs='+', help='CoLoRe files to load (list)')
    parser.add_argument('--valid_pixels', default=None, nargs='+', type=int, help='Downsample in pixels')
    parser.add_argument('--downsampling', default=1, type=float, help='Downsampling to use when extracting objects from the output')
    parser.add_argument('--rsd', action='store_true', help='Apply rsd to the objects')
    parser.add_argument('--simplified', action='store_true', help='Write simplified cats (only RA DEC Z)')
    parser.add_argument('--log-level', default='WARNING', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'NOTSET'])

    args = parser.parse_args()
    level = logging.getLevelName(args.log_level)
    logging.basicConfig(stream=sys.stdout, level=level)

    colore_to_drq_cat(args.in_sim, args.out_file, ifiles=args.ifiles, source=args.source, downsampling=args.downsampling, rsd=args.rsd, valid_pixels=args.valid_pixels, simplified=args.simplified)

def trim_catalog_into_pixels_script(): #pragma: no cover
    parser = argparse.ArgumentParser(description='Method to trim a catalog into multiple catalogs (one for each healpix pixel')

    parser.add_argument('--in-cat', required=True, type=Path, help='Input catalog')
    parser.add_argument('--out-path', required=True, type=Path, help='Output dir for cats (should not exist')
    parser.add_argument('--nside', required=True, type=int, help='Pixelization of the sky (npix = 12*nside**2)')
    parser.add_argument('--log-level', default='WARNING', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'NOTSET'])

    args = parser.parse_args()
    define_logger(args.log_level)

    trim_catalog_into_pixels(args.in_cat, args.out_path, args.nside)

def generate_randoms_from_drq(): #pragma: no cover
    parser = argparse.ArgumentParser(description='Generate a random catalog using an existing one to copy footprint and redshift distribution')

    parser.add_argument('--in-cat', required=True, type=Path, help='Input catalog')
    parser.add_argument('--out-cat', required=True, type=Path, help='Output catalog')
    parser.add_argument('--factor', required=False, type=float, default=1, help='Size factor for the random catalog (compared to the input catalog). (factor=3 means a random catalog three times larger than the input catalog.)')
    parser.add_argument('--pixel_footprint', type=int, nargs='+', default=None, help='Footprint randoms should span. (Default: None, full_sky)')
    parser.add_argument('--nside', required=False, type=int, default=2, help='Pixelizaton of the footprint. (Default 2)')
    parser.add_argument('--method', default='full_sky_dice', choices=['full_sky_dice', 'pixel'], help='Method to generate random positons in the sky: full_sky_dice (throw randoms at the whole sky and then removing the ones not matching) ; pixel (make random pixel by pixel)')
    parser.add_argument('--log-level', default='WARNING', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'NOTSET'])

    args = parser.parse_args()
    define_logger(args.log_level)

    from LyaPlotter.file_types import FilesBase
    in_cat = FilesBase(args.in_cat)

    generate_random_objects_from_data(in_cat.z, in_cat.RA, in_cat.DEC, args.out_cat, args.factor, args.pixel_footprint, args.nside, method=args.method)

def generate_randoms_from_dndz(): #pragma: no cover
    parser = argparse.ArgumentParser(description=' Generate random objects from a given dndz distribution.')

    parser.add_argument('--dndz_file', required=True, type=Path, help='Input filename to read dndz')
    parser.add_argument('--out-cat', required=True, type=Path, help='Output catalog')
    parser.add_argument('--zmin', type=float, default=None, help='Min. redshift. (Default: min redshift in file')
    parser.add_argument('--zmax', type=float, default=None, help='Max. redshift. (Default: max redshift in file')
    parser.add_argument('--factor', required=False, type=float, default=1, help='Size factor for the random catalog (compared to the input catalog). (factor=3 means a random catalog three times larger than the input catalog.)')
    parser.add_argument('--pixel_footprint', type=int, nargs='+', default=None, help='Footprint randoms should span. (Default: None, full_sky)')
    parser.add_argument('--nside', required=False, type=int, default=2, help='Pixelizaton of the footprint. (Default 2)')
    parser.add_argument('--method', default='full_sky_dice', choices=['full_sky_dice', 'pixel'], help='Method to generate random positons in the sky: full_sky_dice (throw randoms at the whole sky and then removing the ones not matching) ; pixel (make random pixel by pixel)')
    parser.add_argument('--log-level', default='WARNING', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'NOTSET'])

    args = parser.parse_args()
    define_logger(args.log_level)

    from LyaPlotter.file_types import FilesBase

    generate_random_objects_from_nz_file(dndz_file=args.dndz_file, out_file=args.out_cat, zmin=args.zmin, zmax=args.zmax, N_factor=args.factor, pixel_footprint=args.pixel_footprint, nside=args.nside, method=args.method)


def define_logger(level):
    level = logging.getLevelName(level)
    logging.basicConfig(stream=sys.stdout, level=level, format='%(levelname)s:%(name)s:%(funcName)s:%(message)s')
