#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
HW3

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
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
import warnings


# CS Ignore warnings
warnings.filterwarnings('ignore')


# CS Use argparse to make an informative help message
parser = argparse.ArgumentParser()
parser.add_argument("Description:", help = "This function examines the u-band properties of Legacy Surveys point sources that have infared excess and are detected in the radio by the FIRST survey and are around (ra,dec) = (163,50)")

# CS Read in FIRST radio data
Fdata = Table.read("/d/scratch/ASTR5160/data/first/first_08jul16.fits")


# CS Assemble a list of all the possible names of the sweep files with glob
fns = glob("/d/scratch/ASTR5160/data/legacysurvey/dr9/north/sweep/9.0/*fits")


boxobj = Table(names=("RA", "DEC")) # Make a box that encloses 3 deg in all directions
boxobj.add_row((166,53))
boxobj.add_row((166,47))
boxobj.add_row((160,47))
boxobj.add_row((160,53))

# CS Determine which FIRST sources lie in the astronomers survey
# CS A circular region centered at ra 163, and dec 50:
ra = 163.
dec = 50.
c_center = SkyCoord(ra*u.degree,dec*u.degree)

# CS Use seperation to get F_coords within 3 deg of 163, 50
fra = np.array(Fdata["RA"])*u.degree
fdec = np.array(Fdata["DEC"])*u.degree
F_coords = SkyCoord(fra,fdec,frame='icrs')
ii = (c_center.separation(F_coords) < 3*u.degree)
sweep_close = Fdata[ii]


# CS Loop through all of the sweep files and use seperation
# CS to get the one that contains our object
sweeplist = [] # CS Global sweep list variable

for i in range(len(fns)):
    ramin, ramax, decmin, decmax = decode_sweep_name(fns[i])
    radecbox = [ramin, ramax, decmin, decmax]
    ii = is_in_box(sweep_close,radecbox) # CS Find the sweep files we need to get all of points in the sweeps that are within 3deg of the location
    if np.any(ii):
        #print(fns[i])
        sweeplist.append(str(fns[i]))


# CS Get every source with r mag < 22 and W1-W2 > 0.5
# CS and is 3 deg of the

rl22 = float(10**(-(22-22.5)/2.5)) # CS Corresponds to flux =  r band mag of 22


# CS Make a Skycoord object that is at the center of 
# CS our circular region, and make it have as many
# CS coords as we have sweeplists

sweep_out = [] # CS Global output variable

# CS Run through each sweep, and apply mag cut
for i in range(len(sweeplist)):
    sweep = Table.read(sweeplist[i], memmap = True) # CS Read in the sweep with memmap=True for speed
    
    # CS Filter bad fluxes
    good_flux = (sweep["FLUX_R"] > 0) & (sweep["FLUX_W1"] > 0) & (sweep["FLUX_W2"] > 0)
    
    # CS Get fluxes
    flux_R = sweep["FLUX_R"]
    flux_w1 = sweep["FLUX_W1"]
    flux_w2 = sweep["FLUX_W2"]
    
    # CS Conver flux to mags (recalling nanomaggies conversions)
    w1_mag = 22.5 - 2.5 * np.log10(sweep["FLUX_W1"]) # CS Corresponds to flux = W1
    w2_mag = 22.5 - 2.5 * np.log10(sweep["FLUX_W2"]) # CS Corresponds to flux = W2
    mag_mask = (sweep["FLUX_R"] > rl22) & ((w1_mag-w2_mag) > 0.5)
    
    sweep = sweep[mag_mask]
    # CS Now make sweep skycoord
    ra_sweep = np.array(sweep["RA"])*u.degree
    dec_sweep = np.array(sweep["DEC"])*u.degree
    c2 = SkyCoord(ra_sweep,dec_sweep, frame='icrs')
    
    # CS Make a mask that corresponds to a flux BRIGHTER than r = 22 and a W1_W2 flux DIMMER than 0.5
    mask_sep = (c_center.separation(c2) < 3*u.degree)
    sweep_out.append(sweep[mask_sep])


# CS Coordinate match FIRST sources and Lecagy Sweep
magcut = vstack(sweep_out)


# CS Employ SkyCoord and search_around_sky at a 1" matching radius
swclra = np.array(sweep_close["RA"])*u.degree
swcldec = np.array(sweep_close["DEC"])*u.degree

c1 = SkyCoord(swclra,swcldec, frame='icrs')


pra2 = np.array(magcut["RA"])*u.degree
pdec2 = np.array(magcut["DEC"])*u.degree
c2 = SkyCoord(pra2,pdec2,frame='icrs')

id1, id2, d2, d3 = c2.search_around_sky(c1, (1/3600)*u.deg) # CS 1" matching radius

first_crossmatch = np.unique(magcut[id1]) # CS Make sure and remove duplicates


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
    # CS Usse subprocess to return output of python sdssDR9query ra dec
    result = subprocess.check_output(["python", "sdssDR9query.py", str(ra), str(dec)],text=True,stderr=subprocess.STDOUT)

    # CS Split output fluxes by comma
    split_result = result.strip().split(",")
    return split_result
   
    
# CS Make global variables to extract ubrite from
u_sdss_crossmatch = []
i_sdss_crossmatch = []
ra_sdss_crossmatch = []
dec_sdss_crossmatch = []

# CS List of coords
coords = list(zip(first_crossmatch["RA"],first_crossmatch["DEC"]))

# CS Use ThreadPoolExecutor to speed up querying
with ThreadPoolExecutor(max_workers=8) as executor:
    # CS Executor.submit schedules a run_sdss_query
    futures = [executor.submit(run_sdss_query, ra, dec) for ra, dec in coords]
    # CS Futures represents the future output from run_sdss_query
    
    for future in as_completed(futures): # CS as_completed only gives us query outputs when they are done, so the code isnt "waiting around"
        sf1 = future.result() # CS Return result of queries
        is_letter(sf1[2], u_sdss_crossmatch) # CS Use is_letter to append only numbers to our tables
        is_letter(sf1[5], i_sdss_crossmatch) # CS This ensures we don't have any not found ones
        is_letter(sf1[0], ra_sdss_crossmatch) # CS Same with ra and dec
        is_letter(sf1[1], dec_sdss_crossmatch)


# CS Get ubrite1 and append
ii = np.argmin(np.unique(u_sdss_crossmatch)) # CS u-band brightest (smallest mag), checking for duplicate sources
ubrite1_u = u_sdss_crossmatch[ii] # CS Get rest of ubrite1 params
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
    #plt.show()


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
ubrite1_u:
19.84357
ubrite1_i:
19.50188
ubrite1_ra:
159.417360038982
ubrite1_dec:
50.465377060102
#3: Number of cross matched r<22 W1-W2 > 0.5 sources:
42
#5 Total number of sources with a match in SDSS:
38
% of #3:
90.47619047619048 %


"""
    
