#!/usr/bin/env python

'''
    Script to obtain the neighbour pixels for a given pixel.
'''
from LyaPlotter.pixel_tools import pixel_get_all_neighbours
import argparse

def main(): #pragma: no cover
    parser = argparse.ArgumentParser(description='Obtain neighbour pixels for a given pixel.')
    parser.add_argument('-n','--nside', required=True, type=int, help="nside for the pixelization of the sky")
    parser.add_argument('-p','--pixel', required=True, type=int, help='pixel from which to obtian neighbours (pixel number')
    parser.add_argument('-d','--depth', required=False, type=int, default=1, help='max depth to get neighbours (e.g: depth=2 will obtain neigbhours of neighbours).')
    parser.add_argument('-r','--rings', required=False, type=bool, default=False, help='whether to output neighbours distributed in rings')

    args = parser.parse_args()

    print(pixel_get_all_neighbours(args.nside, args.pixel, depth=args.depth, rings=args.rings))
    return

if __name__ == '__main__':
    main()
