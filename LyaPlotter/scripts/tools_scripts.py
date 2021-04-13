import argparse
from pathlib import Path
from LyaPlotter.tools import master_to_qso_cat, colore_to_drq_cat
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
    parser.add_argument('--log-level', default='WARNING', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'NOTSET'])

    args = parser.parse_args()
    level = logging.getLevelName(args.log_level)
    logging.basicConfig(stream=sys.stdout, level=level)

    colore_to_drq_cat(args.in_sim, args.out_file, ifiles=args.ifiles, source=args.source, downsampling=args.downsampling, rsd=args.rsd, valid_pixels=args.valid_pixels)