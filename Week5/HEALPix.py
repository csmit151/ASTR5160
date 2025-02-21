#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Areas on the SPhere and HEALPix

Chase L. Smith
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy import units as u
import healpy as hp
from numpy.random import random

#Generate random coordinates
ra = 360.*(random(10**6))
dec = (180/np.pi)*np.arcsin(1.-2*random(10**6))
nside = 1
hpcinside = []


#Use ang2pix with nside=1 to determine which pixels each of your ra and de points
#lie within nside = 1

if __name__ == "__main__":
    hp6 = hp.ang2pix(nside, ra, dec, lonlat=True)
    for i in range(len(hp6)):
        if hp6[i] == 0:
            hpcinside.append(hp6[i])
    print("Number of pixels within nside=1",len(hpcinside))
    print("The area of a nside=1 HEALPixel is:",hp.nside2pixarea(2**0),"square radians")
    print("Number of points in each HEAlPixel:",np.unique(hp6,return_counts=True))
    
    """
    Number of pixels within nside=1 83093
    The area of a nside=1 HEALPixel is: 1.0471975511965976 square radians
    Number of points in each HEAlPixel: (array([ 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11]), array([83093, 83772, 83309, 83132, 83319, 83627, 83385, 83206, 83541,
       83239, 83145, 83232]))
    So, the results are consistent with the number of pixels being equal in area
    """
    

