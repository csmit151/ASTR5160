#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
HW4

@author: csmit151
"""

# CS Import relevant packages
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack
import os
from tqdm import tqdm
from Week8.CrossMatch import decode_sweep_name
from Week8 import sdssDR9query
from Week8.CrossMatch import is_in_box
import argparse
from glob import glob
from astropy.coordinates import SkyCoord
from astropy import units as u
import warnings


# CS Ignore warnings
warnings.filterwarnings('ignore')

# CS Step 1 is to determine cuts based off of
# CS the information in qs316

"""
MATCHING SECTION
"""

# CS Read in 316 confirmed quasars
qs316 = Table.read("/d/scratch/ASTR5160/week10/qsos-ra180-dec30-rad3.fits", memmap=True)

#print(qs316.colnames)
# CS ['NAME', 'RA', 'DEC', 'ZEM', 'GMAG', 'SOURCE']
# CS ZEM = emmission redshift
# CS Source is as: SDSSDR7QSO

# CS Let's start by crossmatching these quasars to 
# CS sources in the sweep files to get more informatin
# CS about over-arcing trends present

# CS Assemble a list of all the possible names of the sweep files with glob
fns = glob("/d/scratch/ASTR5160/data/legacysurvey/dr9/south/sweep/9.0/*fits")

# CS A circular region centered at ra 180, dec 30:
ra = 180.
dec = 30.
c_center = SkyCoord(ra*u.degree,dec*u.degree)
rl19 = float(10**(-(19-22.5)/2.5)) # CS Corresponds to flux =  r band mag of 22


sweeplist_in = [] # CS Global sweep list variable

for i in range(len(fns)):
    ramin, ramax, decmin, decmax = decode_sweep_name(fns[i])
    radecbox = [ramin, ramax, decmin, decmax]
    ii = is_in_box(qs316,radecbox) # CS Find the sweep files we need to get all of points in the sweeps 
    if np.any(ii):
        sweeplist_in.append(str(fns[i]))

print("Sweep List in:",sweeplist_in)
print("len(sweeplist_in):",len(sweeplist_in))

sweep_tables = []
for i in tqdm(range(len(sweeplist_in))):
    # CS Put sweep file objects in one table
    sweep = Table.read(sweeplist_in[i], memmap=True)
    mag_mask = (sweep["FLUX_R"] > rl19)
    sweep = sweep[mag_mask]
    # CS Now make sweep skycoord
    ra_sweep = np.array(sweep["RA"])*u.degree
    dec_sweep = np.array(sweep["DEC"])*u.degree
    c2 = SkyCoord(ra_sweep,dec_sweep, frame='icrs')
    # CS Use separation
    mask_sep = (c_center.separation(c2) < 3*u.degree)
    sweep_tables.append(sweep)

# CS Combine sweep files
sw_stack = vstack(sweep_tables)

# CS Coordinate match sweep sources with our confirmed quasars
# CS Employ SkyCoord and search_around_sky at a 1" matching radius
qsra = np.array(qs316["RA"])*u.degree # CS 316 Quasars
qsdec = np.array(qs316["DEC"])*u.degree
c1 = SkyCoord(qsra,qsdec, frame='icrs')

sweep_ra = np.array(sw_stack["RA"])*u.degree # CS our combined sweep files
sweep_dec = np.array(sw_stack["DEC"])*u.degree # CS wtih coords < 3 deg of center
c2 = SkyCoord(sweep_ra,sweep_dec,frame='icrs')
id1, id2, d2, d3 = c1.search_around_sky(c2, (1/3600)*u.deg) # CS 1" matching radius

# CS These are backwards but code crashes without this config?
qs316_match = qs316[id2] # CS quasar files matched objects
qs_crossmatch = sw_stack[id1] # CS Sweep files matched objects


print("qs_crossmatch len:",len(qs_crossmatch["RA"]))
print("qs316_match len:",len(qs316_match))
print("max(qs_crossmatch[RA]):",np.max(qs_crossmatch["RA"]))
# CS max(qs_crossmatch[RA]) is too large....

# CS Hmm prints out 190....

# CS So now qs_crossmatch is just quasars, but we have all
# CS of the info on them
# CS Now lets look for trends

#print(qs316_match.colnames)
#print(qs_crossmatch.colnames)
"""
g_z = []
r_w1 = []

fig, axs = plt.subplots(1, 1, figsize=(8,6))

for i in range(len(qs_crossmatch["FLUX_G"])):
    axs.scatter(i,qs_crossmatch["FLUX_G"][i]/qs_crossmatch["MW_TRANSMISSION_G"][i])
    axs.scatter(i,qs_crossmatch["FLUX_R"][i]/qs_crossmatch["MW_TRANSMISSION_R"][i])
    axs.scatter(i,qs_crossmatch["FLUX_Z"][i]/qs_crossmatch["MW_TRANSMISSION_Z"][i])
    axs.scatter(i,qs_crossmatch["FLUX_W1"][i]/qs_crossmatch["MW_TRANSMISSION_W1"][i])
    axs.scatter(i,qs_crossmatch["FLUX_W2"][i]/qs_crossmatch["MW_TRANSMISSION_W2"][i])
    axs.scatter(i,qs_crossmatch["FLUX_W3"][i]/qs_crossmatch["MW_TRANSMISSION_W3"][i])
    axs.scatter(i,qs_crossmatch["FLUX_W4"][i]/qs_crossmatch["MW_TRANSMISSION_W4"][i])
    #plt.scatter(i,qs316_match["ZEM"][i])
    axs.set_xlabel("Object")
    axs.set_ylabel("Flux")
    
plt.savefig("FluxTest.png")
"""

"""
FUNCTION SECTION
"""

def splendid_function(T_in):
    print("My table is:",T_in)
    
    
# CS args.Table = table input as python HW4.py mTable.fits
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("Table")
    args = parser.parse_args()
    splendid_function(args.Table)

