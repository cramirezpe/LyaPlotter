'''
    Module with tools related to healpy pixels.
'''

import healpy as hp
import numpy as np

def pixel_get_all_neighbours(nside, pixel, depth=1, rings=False, nest=True):
    ''' 
        Function to get all neighbour pixels to a given depth. (Included the given one)

        Args: 
            nside (int): pixelization of the sky.
            pixel (int): pixel number from which obtain neighbours.
            depth (int,optional): max depth to get neighbours (e.g: depth=2 will obtain neighbours of neighbours). (default: 1)
            rings (bool, optional): whether to output neighbours distributed in rings (default: False)
            nest (bool, optional): if True, assume NESTED pixel ordering (the one used in LyaCoLoRe by default), otherwise, RING pixel ordering.

        Returns:
            set of neighbour pixels.
    '''

    all_elements = set([pixel])
    neighbours = []

    neighbours.append(set([pixel]))
    all_elements.update(neighbours[0])

    for i in range(depth):
        neighbours.append(set())
        for pix in neighbours[i]:
            new_elements = set(hp.pixelfunc.get_all_neighbours(nside,pix, nest=nest))
            new_elements.difference_update(all_elements)
            neighbours[i+1].update(new_elements)
            all_elements.update(new_elements)
        
    if rings:
        return neighbours
    else:
        return all_elements

def plot_pixels(nside, pixels, nest=True, ax=False): #pragma: no cover
    ''' 
        Plot the given pixels in a mollview. 
    '''
    if nest:
        pixels = hp.pixelfunc.nest2ring(16,list(pixels))

    m = np.full(12*nside**2, 0, dtype=np.double)
    for i in pixels:
        n = i -1
        m[n] = 1

    if ax: 
        plt.axes(ax)
        hold = True
    else:
        hold = False
    hp.mollview(m, hold=hold)

