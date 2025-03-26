#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Classification in Imaging

@author: csmit151
"""

import numpy as np
import matplotlib.pyplot as plt
from Week8.CrossMatch import is_in_box
from Week8.CrossMatch import decode_sweep_name
from astropy.table import Table, join
from tqdm import tqdm


stars = Table.read("/d/scratch/ASTR5160/week10/stars-ra180-dec30-rad3.fits")
qs = Table.read("/d/scratch/ASTR5160/week10/qsos-ra180-dec30-rad3.fits")

def whichSweepSouth(objects):
    """ Takes an array of South RA/Dec locations and returns
    which sweep files would need to be read to find objects
    at those locations.
    
    Parameters
    ----------
    ra: :class: 'array'
        An array of RA values
    dec: :class: 'array'
        An array of DEC values
        
    
    Returns
    -------
    List of sweep files
    """
    # CS start by grabbing the list of sweep boxes we want to search in:
    sweepradecs = []
    # CS iterating through all avaliable sweep files
    # CS and getting the ra and dec values from them
    for sweepname in os.listdir("/d/scratch/ASTR5160/data/legacysurvey/dr9/south/sweep/9.0"):
        filepath = os.path.join("/d/scratch/ASTR5160/data/legacysurvey/dr9/south/sweep/9.0", sweepname)
        if os.path.isfile(filepath):
            if filepath != "/d/scratch/ASTR5160/data/legacysurvey/dr9/south/sweep/9.0/legacysurvey_dr9_south_sweep_9.0.sha256sum":
                sweepradecs.append(decode_sweep_name(sweepname))
    # CS now determine whether or not a given object is in a given box
    sweeps = []
    swradecs = []
    for i in range(len(sweepradecs)):
        for k in range(len(objects)):
            sweeps.append(is_in_box(objects[k],sweepradecs[i]))
            swradecs.append(sweepradecs[i])
    sweepscor = []
    swradecscor = []
    
    
    for i in range(len(sweeps)):
        if sweeps[i] == True:
            sweepscor.append(sweeps[i])
            swradecscor.append(swradecs[i])
    print("I recover:",len(np.unique(swradecscor, axis=0)),"sweep files")
    print(np.unique(swradecscor, axis=0))

"""
Coord matched star sweep files:
('I recover:', 4, 'sweep files')
[[170. 180.  25.  30.]
 [170. 180.  30.  35.]
 [180. 190.  25.  30.]
 [180. 190.  30.  35.]]
Coord matched qs sweep files:
('I recover:', 4, 'sweep files')
[[170. 180.  25.  30.]
 [170. 180.  30.  35.]
 [180. 190.  25.  30.]
 [180. 190.  30.  35.]]
"""    



if __name__ == "__main__":
    print("Coord matched star sweep files:")
    whichSweepSouth(stars)
    print("Coord matched qs sweep files:")
    whichSweepSouth(qs)