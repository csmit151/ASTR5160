#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Flagging Bad Data in Imaging

@author: csmit151
"""


# CS Import relevant packages
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, join
import os
from tqdm import tqdm
from glob import glob
from Week8.CrossMatch import decode_sweep_name
from Week8.CrossMatch import is_in_box
from astropy.coordinates import SkyCoord
from astropy import units as u


"""
flag = 2**2 # CS Corresponds to saturation in g band
ii = (objs["MASKBITS"] & flag) == 0
goodobjs = objs[ii]
"""

# CS Find the closest object in the sweep files to 
ra = 188.53667
dec = 21.04572

objs = Table(names=("RA", "DEC"))
objs.add_row((ra,dec))

fns = glob("/d/scratch/ASTR5160/data/legacysurvey/dr9/south/sweep/9.0/*fits")

qs_to_match = Table.read("/d/scratch/ASTR5160/week10/qsos-ra180-dec30-rad3.fits")


"""
# CS Loop through all of the sweep files and use is_in_box to
# CS  get the one that contains our object6
for i in range(len(fns)):
    ramin, ramax, decmin, decmax = decode_sweep_name(fns[i])
    radecbox = [ramin, ramax, decmin, decmax]
    ii = is_in_box(objs,radecbox)
    if np.any(ii):
        print(fns[i])
""" 

"""
CS Our sweep file is:
    /d/scratch/ASTR5160/data/legacysurvey/dr9/south/sweep/9.0/sweep-180p020-190p025.fits
"""

#sweep1 = Table.read("/d/scratch/ASTR5160/data/legacysurvey/dr9/south/sweep/9.0/sweep-180p020-190p025.fits")

#ii = (sweep1["RA"] > (ra-0.0001)) & (sweep1["RA"] < (ra+0.0001)) & (sweep1["DEC"] > (dec-0.0001)) & (sweep1["DEC"] < (dec+0.0001))

#smallsweep = sweep1[ii]


# CS Retrive every point source in the sweeps that are within 3 deg of
# CS ra, dec = 180, 30

# CS Star by finding which sweep files we need
# CS Loop through all of the sweep files and use is_in_box to
# CS  get the one is in our box

boxobj = Table(names=("RA", "DEC")) # Make a box that encloses 3 deg in all directions
boxobj.add_row((183,33))
boxobj.add_row((177,33))
boxobj.add_row((177,27))
boxobj.add_row((183,27))




sweeplist = []

for i in range(len(fns)):
    ramin, ramax, decmin, decmax = decode_sweep_name(fns[i])
    radecbox = [ramin, ramax, decmin, decmax]
    ii = is_in_box(boxobj,radecbox) # CS Find the sweep files we need to get all of points in the sweeps that are within 3deg of the location
    if np.any(ii):
        #print(fns[i])
        sweeplist.append(str(fns[i]))

# CS Get every point source with TYPE == PSF and r mag < 20
sw1 = []
sw2 = []
sw3 = []
sw4 = []

swtotlen = 0

for i in tqdm(range(len(sweeplist))):
    ra1, dec1 = [180.]*u.degree,[30.]*u.degree
    c1 = SkyCoord(ra1,dec1,frame='icrs')
    if i == 0:
        sweep = Table.read(sweeplist[i])
        ra2 = np.array(sweep["RA"])*u.degree
        dec2 = np.array(sweep["DEC"])*u.degree
        c2 = SkyCoord(ra2,dec2,frame='icrs')
        kk = (sweep["FLUX_R"] > 10) & (c1.separation(c2) < 3*u.degree)# CS Corresponds to r < 20
        swa = sweep[kk]
        ii = (swa["TYPE"] == "PSF")
        sw1 = swa[ii]
        print(len(sw1["RA"]))
        swtotlen = swtotlen + len(sw1["RA"])
    if i == 1:
        sweep = Table.read(sweeplist[i])
        ra2 = np.array(sweep["RA"])*u.degree
        dec2 = np.array(sweep["DEC"])*u.degree
        c2 = SkyCoord(ra2,dec2,frame='icrs')
        kk = (sweep["FLUX_R"] > 10) & (c1.separation(c2) < 3*u.degree)# CS Corresponds to r < 20
        swa = sweep[kk]
        ii = (swa["TYPE"] == "PSF")
        sw2 = swa[ii]
        print(len(sw2["RA"]))
        swtotlen = swtotlen + len(sw2["RA"])
    if i == 2:
        sweep = Table.read(sweeplist[i])
        ra2 = np.array(sweep["RA"])*u.degree
        dec2 = np.array(sweep["DEC"])*u.degree
        c2 = SkyCoord(ra2,dec2,frame='icrs')
        kk = (sweep["FLUX_R"] > 10) & (c1.separation(c2) < 3*u.degree)# CS Corresponds to r < 20
        swa = sweep[kk]
        ii = (swa["TYPE"] == "PSF")
        sw3 = swa[ii]
        print(len(sw3["RA"]))
        swtotlen = swtotlen + len(sw3["RA"])
    if i == 3:
        sweep = Table.read(sweeplist[i])
        ra2 = np.array(sweep["RA"])*u.degree
        dec2 = np.array(sweep["DEC"])*u.degree
        c2 = SkyCoord(ra2,dec2,frame='icrs')
        kk = (sweep["FLUX_R"] > 10) & (c1.separation(c2) < 3*u.degree)# CS Corresponds to r < 20
        swa = sweep[kk]
        ii = (swa["TYPE"] == "PSF")
        sw4 = swa[ii]
        print(len(sw4["RA"]))
        swtotlen = swtotlen + len(sw4["RA"])


print(swtotlen)

# CS Coordinate match psfobjs to the qsosfile

psfobjs = np.concatenate((sw1,sw2,sw3,sw4), axis=0)

qra1 = np.array(qs_to_match["RA"])*u.degree
qdec1 = np.array(qs_to_match["DEC"])*u.degree
c1 = SkyCoord(qra1,qdec1)

pra2 = np.array(psfobjs["RA"])*u.degree
pdec2 = np.array(psfobjs["DEC"])*u.degree
c2 = SkyCoord(pra2,pdec2,frame='icrs')

id1, id2, d2, d3 = c2.search_around_sky(c1, (1/3600)*u.deg) # CS 1" mr

# CS Bah ha! It works



if __name__ == "__main__":
    #print(len(psfobjs))
    #print(swtotlen)
    print("Number of cross matched r<20 qsos:")
    print(len(qra1[id1]))
    

    """
    print(smallsweep["RA"][0])
    print(smallsweep["DEC"][0])
    print(smallsweep["TYPE"][0])
    print(smallsweep["OBJID"][0])
    print(smallsweep["FLUX_G"][0])
    print(smallsweep["FLUX_R"][0])
    print(smallsweep["FLUX_Z"][0])

    flag = 2**1 # CS Corresponds to saturation in g band
    ii = (smallsweep["ALLMASK_G"][0] & flag) == 0
    print(ii)
    if np.any(ii):
        print("Not saturated in G")
    else:
        print("Saturated in G")
        
    flag = 2**1 # CS Corresponds to saturation in r band
    jj = (smallsweep["ALLMASK_R"][0] & flag) == 0
    print(jj)
    if np.any(jj):
        print("Not saturated in R")
    else:
        print("Saturated in R")
        
    flag = 2**1 # CS Corresponds to saturation in z band
    kk = (smallsweep["ALLMASK_Z"][0] & flag) == 0
    print(kk)
    if np.any(kk):
        print("Not saturated in Z")
    else:
        print("Saturated in Z")
    """


"""
RESULTS
-------------------

RA: 188.53664112612532
DEC: 21.045714562398302
Obj Type: PSF
Obj ID: 4782
G:1155.0397
R:1747.2861
Z:0.0
Saturated in G
Saturated in R
Saturated in Z
This object is saturated in all three bands
This object is saturated in Legacy Surveys Sky Viewer, and is not a galxay

Total number of PSF objects: 39655
Number of cross matched r<20 qsos:
275
"""
