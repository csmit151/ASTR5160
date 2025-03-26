#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Magnitude Systems

@author: csmit151
"""

import numpy as np
import matplotlib.pyplot as plt
from Week8.CrossMatch import is_in_box
from Week8.CrossMatch import decode_sweep_name
from astropy.table import Table, join
from astropy.io import fits
from tqdm import tqdm
import os
import math

"""
CS For PG1633+099A, convert UBVRI into ugriz,
and show g and z mags.

DR9 g = 15.70
DR9 z = 14.55
"""


"""
# CS Calculate UBVRI
V = 15.256
B_V = 0.873
U_B = 0.320
V_R = 0.505
R_I = 0.511
B = B_V + V
U = U_B + B
I = -(V_R+R_I)+V
R = R_I+I

# CS Utilize Jester 2005 R_I < 1.15:
g = V +0.60*(B_V) - 0.12
r_z = 1.72*(R_I) - 0.41
r = V - 0.46*(B_V)+0.11
z = -(r_z)+r

if __name__ == "__main__":
    print("For V:",V,", B_V:", B_V, ", and R_I:", R_I, "I find g:", g, " and z:", z)
    print("DR9 g = 15.70, DR9 z = 14.55")
"""

"""
CS For V: 15.256 , B_V: 0.873 , and R_I: 0.511 I find g: 15.6598  and z: 14.4955, DR9 g = 15.70, DR9 z = 14.55
So pretty darn close!
""" 

# CS Obtain the fluxes for PG1633+099A from the Legacy Surveys sweep files

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
    print(swradecscor)

# CS coords of PG1633+099A:
# 248.8583333    9.7980556

# CS so is in sweep-240p005-250p010.fits
# J163525.98+094753.1

sweep = "/d/scratch/ASTR5160/data/legacysurvey/dr9/south/sweep/9.0/sweep-240p005-250p010.fits"

with fits.open(sweep) as hdul:
    data = hdul[1].data
    mask = (data["RA"] > 248.8582) & (data["RA"] < 248.8583) & (data["DEC"] > 9.79) & (data["DEC"] < 9.80)
    PGA = data[mask]
    #print(data.dtype.names)
    print("RA:",PGA["RA"])
    print("DEC:",PGA["DEC"])
    print("PAR:",PGA["PARALLAX"])
    print("FG:",PGA["FLUX_G"])
    print("FR:",PGA["FLUX_R"])
    print("FZ:",PGA["FLUX_Z"])
    print("W1:",PGA["FLUX_W1"])
    print("W2:",PGA["FLUX_W2"])
    print("W3:",PGA["FLUX_W3"])
    print("W4:",PGA["FLUX_W4"])


# CS Navigate Values
g = 15.70
r = 15.19
z = 14.55

# CS Fluxes from sweep file (in nanomaggies)
FG = 586.09576
FR = 1094.6772
FZ = 1559.2128

print("g mag:",22.5-2.5*math.log10(FG))
print("r mag:",22.5-2.5*math.log10(FR))
print("z mag:",22.5-2.5*math.log10(FZ))

# CS Wise band mags:
print("W1 mag:",22.5-2.5*math.log10(603.8309))
print("W2 mag:",22.5-2.5*math.log10(328.57047))
print("W3 mag:",22.5-2.5*math.log10(38.98512))
print("W4 mag: No Detection (Negative Flux)")

"""
('g mag:', 15.580078551072756)
('r mag:', 14.90178481833016)
('z mag:', 14.517736521587057)
CS Pretty close to navigate tool values!
('W1 mag:', 15.54771166619171)
('W2 mag:', 16.208428677513645)
('W3 mag:', 18.522752811601237)
W4 mag: No Detection (Negative Flux)
"""

