#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
CrossMatching Surveys

"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
import os
from tqdm import tqdm


# CS read in and plot first data
objsVLA = Table.read("/d/scratch/ASTR5160/data/first/first_08jul16.fits")
#plt.scatter(objsVLA["RA"],objsVLA["DEC"], marker = "o")




# CS Write Python code that reads in the FIRST radio data
# CS takes the coords of the first 100 objs and uses
# CS sdssDR9query.py to obtain the SDSS optical data

# CS first 100 objects in FIRST
f100ra = str(objsVLA["RA"][:100])
f100dec = str(objsVLA["DEC"][:100])


"""
print(str(objsVLA["RA"][0]))

for i in tqdm(range(len(f100ra))):
    os.system("python sdssDR9query.py " + str(objsVLA["RA"][i]) +" " + str(objsVLA["DEC"][i]) + " >> F100.txt")
# CS very very very slow


# CS read in and plot more recent data
#objsVLA = Table.read("/d/scratch/ASTR5160/data/legacysurvey/dr9/north/sweep/9.0/sweep-320m005-330p000.fits")
#plt.scatter(objsVLA["RA"],objsVLA["DEC"], marker = "o")
"""

# CS imorted from WyoCourses via ADM

def is_in_box(objs, radecbox):
    """Determine which of an array of objects are inside an RA, Dec box.

    Parameters
    ----------
    objs : :class:`~numpy.ndarray`
        An array of objects. Must include at least the columns "RA" and "DEC".
    radecbox : :class:`list`
        4-entry list of coordinates [ramin, ramax, decmin, decmax] forming the
        edges of a box in RA/Dec (degrees).

    Returns
    -------
    :class:`~numpy.ndarray`
        ``True`` for objects in the box, ``False`` for objects outside of the box.

    Notes
    -----
        - Tests the minimum RA/Dec with >= and the maximum with <
    """
    ramin, ramax, decmin, decmax = radecbox

    # ADM check for some common mistakes.
    if decmin < -90. or decmax > 90. or decmax <= decmin or ramax <= ramin:
        msg = "Strange input: [ramin, ramax, decmin, decmax] = {}".format(radecbox)
        log.critical(msg)
        raise ValueError(msg)

    ii = ((objs["RA"] >= ramin) & (objs["RA"] < ramax)
          & (objs["DEC"] >= decmin) & (objs["DEC"] < decmax))

    return ii

# CS imorted from WyoCourses via ADM
def decode_sweep_name(sweepname, nside=None, inclusive=True, fact=4):
    """Retrieve RA/Dec edges from a full directory path to a sweep file

    Parameters
    ----------
    sweepname : :class:`str`
        Full path to a sweep file, e.g., /a/b/c/sweep-350m005-360p005.fits
    nside : :class:`int`, optional, defaults to None
        (NESTED) HEALPixel nside
    inclusive : :class:`book`, optional, defaults to ``True``
        see documentation for `healpy.query_polygon()`
    fact : :class:`int`, optional defaults to 4
        see documentation for `healpy.query_polygon()`

    Returns
    -------
    :class:`list` (if nside is None)
        A 4-entry list of the edges of the region covered by the sweeps file
        in the form [RAmin, RAmax, DECmin, DECmax]
        For the above example this would be [350., 360., -5., 5.]
    :class:`list` (if nside is not None)
        A list of HEALPixels that touch the  files at the passed `nside`
        For the above example this would be [16, 17, 18, 19]
    """
    # ADM extract just the file part of the name.
    sweepname = os.path.basename(sweepname)

    # ADM the RA/Dec edges.
    ramin, ramax = float(sweepname[6:9]), float(sweepname[14:17])
    decmin, decmax = float(sweepname[10:13]), float(sweepname[18:21])

    # ADM flip the signs on the DECs, if needed.
    if sweepname[9] == 'm':
        decmin *= -1
    if sweepname[17] == 'm':
        decmax *= -1

    if nside is None:
        return [ramin, ramax, decmin, decmax]

    pixnum = hp_in_box(nside, [ramin, ramax, decmin, decmax],
                       inclusive=inclusive, fact=fact)

    return pixnum



def whichSweep(objects):
    """ Takes an array of RA/Dec locations and returns
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
    for sweepname in os.listdir("/d/scratch/ASTR5160/data/legacysurvey/dr9/north/sweep/9.0"):
        filepath = os.path.join("/d/scratch/ASTR5160/data/legacysurvey/dr9/north/sweep/9.0", sweepname)
        if os.path.isfile(filepath):
            if filepath != "/d/scratch/ASTR5160/data/legacysurvey/dr9/north/sweep/9.0/legacysurvey_dr9_north_sweep_9.0.sha256sum":
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
    
    
if __name__ == "__main__":
    whichSweep(objsVLA[:100])

