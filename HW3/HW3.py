#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
HW3

@author: csmit151
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack
import os
from tqdm import tqdm
from Week8.CrossMatch import decode_sweep_name
from Week8 import sdssDR9query
import argparse
from glob import glob
from astropy.coordinates import SkyCoord
from astropy import units as u
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed



# CS Use argparse to make an informative help message
parser = argparse.ArgumentParser()
parser.add_argument("Description:", help = "This function examines the u-band properties of Legacy Surveys point sources that have infared excess and are detected in the radio by the FIRST survey and are around (ra,dec) = (163,50)")

# CS Read in the FIRST sources
first_to_match = Table.read("/d/scratch/ASTR5160/data/first/first_08jul16.fits")

# CS Assemble a list of all the possible names of the sweep files with glob
fns = glob("/d/scratch/ASTR5160/data/legacysurvey/dr9/north/sweep/9.0/*fits")

# CS A circular region centered at ra 163, and dec 50
ra = 163.
dec = 50.

# CS Determine which FIRST sources lie in the astronomers survey
c_center = SkyCoord(ra*u.degree,dec*u.degree)

# CS Read in FIRST radio data
Fdata = Table.read("/d/scratch/ASTR5160/data/first/first_08jul16.fits")

# CS Loop through all of the sweep files and use seperation
# CS to get the one that contains our object
sweeplist = [] # CS Global sweep list variable

for i in range(len(fns)):
    ramin, ramax, decmin, decmax = decode_sweep_name(fns[i])
    # CS Convert the center point of the sweep file we're testing
    # CS to a SkyCoord for .separation
    sweep_skycoord = SkyCoord(((ramin+ramax)/2)*u.degree,((decmin+decmax)/2)*u.degree,frame='icrs')
    if c_center.separation(sweep_skycoord) < 3*u.degree:
        sweeplist.append(fns[i]) # CS Only append if we found one inside of 3 deg

# CS Get every source with r mag < 22 and W1-W2 > 0.5
# CS and is 3 deg of the

rl22 = float(10**(-(22-22.5)/2.5)) # CS Corresponds to flux =  r band mag of 22

w1_w2_p05 = float(10**(-(0.5-22.5)/2.5)) # CS Corresponds to flux = W1_W2 of 0.5


# CS Make a Skycoord object that is at the center of 
# CS our circular region, and make it have as many
# CS coords as we have sweeplists

c1 = SkyCoord([ra]*len(sweeplist)*u.degree,[dec]*len(sweeplist)*u.degree,frame='icrs')

sweep_out = [] # CS Global output variable

# CS Run through each sweep, and apply mag cut
for i, sweep_path in enumerate(tqdm(sweeplist)):
    sweep = Table.read(sweep_path) # CS Read in the sweep
    ra_sweep = np.array(sweep["RA"])*u.degree
    dec_sweep = np.array(sweep["DEC"])*u.degree
    c2 = SkyCoord(ra_sweep,dec_sweep, frame='icrs')
    
    sep = c1[i].separation(c2) # CS Use astropy.coordinates seperation from Class 8
    # CS Make a mask that corresponds to a flux BRIGHTER than r = 22 and a W1_W2 flux DIMMER than 0.5
    mask_mag = ((sweep["FLUX_R"] > rl22) & ((sweep["FLUX_W1"] - sweep["FLUX_W2"]) < w1_w2_p05) & (sep < 3*u.degree))
    sweep_out.append(sweep[mask_mag])

# CS Coordinate match FIRST sources and Lecagy Sweep
magcut = vstack(sweep_out)


# CS Employ SkyCoord and search_around_sky at a 1" matching radius
qra1 = np.array(first_to_match["RA"])*u.degree
qdec1 = np.array(first_to_match["DEC"])*u.degree
c1 = SkyCoord(qra1,qdec1)

pra2 = np.array(magcut["RA"])*u.degree
pdec2 = np.array(magcut["DEC"])*u.degree
c2 = SkyCoord(pra2,pdec2,frame='icrs')

id1, id2, d2, d3 = c2.search_around_sky(c1, (1/3600)*u.deg) # CS 1" matching radius

# CS Get sources
first_crossmatch = magcut[id1]

def is_letter(test, list_in):
    """Takes an input and determines
    if it is a number or a letter, and if
    it is a number and can be float() it appends
    it to the list_in
    
    Parameters
    ----------
    test: :class: 'string'
        Unknown letter/number the function is testing
    list_in: :class: 'list'
        Input list that the function will append to
    Returns
    -------
    A list appended with only float() values
    """
    try:
        list_in.append(float(test))
    except ValueError:
        pass


# CS Use sdssDR9query.py to fetch 
# CS SDSS u-band and i-band magnitudes for the cources

def run_sdss_query(ra, dec):
    """ Utilizes subprocess to run
    sdssDR9query.py and split the result by ","
    and then return it
    
    Parameters
    ----------
    ra: :class: 'string'
        ra to input into sdssDR9query.py
    dec: :class: 'string'
        dec to input into sdssDR9query.py
    Returns
    -------
    Resulting split output from sdssDR9query
    """
    result = subprocess.check_output(
        ["python", "sdssDR9query.py", str(ra), str(dec)],
        text=True,
        stderr=subprocess.STDOUT
    )

    # CS Split fluxes by comma
    split_result = result.strip().split(",")
    return split_result
   
    
# CS Make global variables to extract ubrite from
u_sdss_crossmatch = []
i_sdss_crossmatch = []
ra_sdss_crossmatch = []
dec_sdss_crossmatch = []

coords = list(zip(first_crossmatch["RA"],first_crossmatch["DEC"]))

# CS Use ThreadPoolExecutor to speed up querying
with ThreadPoolExecutor(max_workers=8) as executor:
    futures = [executor.submit(run_sdss_query, ra, dec) for ra, dec in coords]

    for future in tqdm(as_completed(futures), total=len(futures)):
        sf1 = future.result()
        is_letter(sf1[2], u_sdss_crossmatch)
        is_letter(sf1[5], i_sdss_crossmatch)
        is_letter(sf1[0], ra_sdss_crossmatch)
        is_letter(sf1[1], dec_sdss_crossmatch)


# CS Get ubrite1 and append
ii = np.argmin(u_sdss_crossmatch) # CS u-band brightest (smallest mag)
ubrite1_u = u_sdss_crossmatch[ii]
ubrite1_i = i_sdss_crossmatch[ii]
ubrite1_ra = ra_sdss_crossmatch[ii]
ubrite1_dec = dec_sdss_crossmatch[ii]


print("ubrite1_u:")
print(ubrite1_u)
print("ubrite1_i:")
print(ubrite1_i)
print("ubrite1_ra:")
print(ubrite1_ra)
print("ubrite1_dec:")
print(ubrite1_dec)

# CS Convert the u-band and i-band mags to fluxes, assuming 
# CS that these fluxes are in nanomaggies


u_flux = 10**((-(float(ubrite1_u)-22.5))/2.5)
i_flux = 10**((-(float(ubrite1_i)-22.5))/2.5)

# CS Plot 9 available fluxes for ubrite1 as a function of wavelength

# CS Employ SkyCoord and search_around_sky at a 1" matching radius
# CS to grab all of the fluxes of ubrite1
ura1 = np.array(ubrite1_ra)*u.degree
udec1 = np.array(ubrite1_dec)*u.degree
c3 = SkyCoord(ura1,udec1)

pra2 = np.array(magcut["RA"])*u.degree
pdec2 = np.array(magcut["DEC"])*u.degree
c2 = SkyCoord(pra2,pdec2,frame='icrs')

qq = c3.separation(c2) < (1/3600)*u.degree # CS 1" matching radius

ubrite1 = magcut[qq] # CS ubrite1 Object

# CS Make an array with the wavelengths of all of our bands in Angstroms
# CS Corresponds to ugrizW1,2,3,4
wavelengths = [3543, 4770, 6231, 7625, 9134, 34000, 46000, 12000, 22000]

# CS Make an array with all of our fluxes
# CS Also corresponds to ugrizW1,2,3,4
fluxes = [u_flux, ubrite1["FLUX_G"][0], ubrite1["FLUX_R"][0],i_flux,ubrite1["FLUX_Z"][0],ubrite1["FLUX_W1"][0],ubrite1["FLUX_W2"][0],ubrite1["FLUX_W3"][0],ubrite1["FLUX_W4"][0]]


# CS Plot wavelength vs flux
def lam_v_flux(lam,flux):
    """ Takes an array of wavelengths and an
    array of fluxes and plots flux as a function
    of wavelength
    
    Parameters
    ----------
    lam: :class: 'array'
        Array of wavelengths
    flux: :class: 'array'
        Array of fluxes
    Returns
    -------
    Plot of flux as a function of wavelength
    """
    plt.scatter(lam,flux)
    plt.xlim(3000,25000)
    plt.xlabel("Band Wavelength (Angstroms)")
    plt.ylabel("Nanomaggies (3631 * $10^9$) Jy")
    plt.show()


def print_results():
    """ This function prints out
    the desired results from the cross matching
    
    Parameters
    ----------
    None
    Returns
    -------
    An annotated list of the HW3 results for easy grading
    """
    print("#3: Number of cross matched r<22 W1-W2 > 0.5 sources:")
    print(len(first_crossmatch["RA"]))
    print("#5 Total number of sources with a match in SDSS:")
    print(len(u_sdss_crossmatch))
    print("% of #3:")
    print(float(len(u_sdss_crossmatch)/len(first_crossmatch["RA"]))*100,"%")
    lam_v_flux(wavelengths,fluxes)


if __name__ == "__main__":
    print_results()
    
"""
Results
----------------
#3: Number of cross matched r<22 W1-W2 > 0.5 sources:
585
#5 Total number of sources with a match in SDSS:
567
So ~96.92% recovered
1brite1 u band: 14.8
ubrite1 i band: 13.63

"""
    
