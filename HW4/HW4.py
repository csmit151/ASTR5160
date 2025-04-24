#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
HW4

@author: csmit151
"""

# CS Import relevant packages
import numpy as np
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
import stat

# CS Ignore warnings
warnings.filterwarnings('ignore')

# CS Step 1 is to determine cuts based off of
# CS the information in qs316

"""
MATCHING SECTION
"""
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
rl19 = float(10**(-(19-22.5)/2.5)) # CS Corresponds to flux =  r band mag of 19
#print("Should be 19:", 22.5-2.5*np.log10(rl19))
#CS It is 19, good to see

sweeplist_in = [] # CS Global sweep list variable

for fns in fns:
    ramin, ramax, decmin, decmax = decode_sweep_name(fns)
    radecbox = [ramin, ramax, decmin, decmax]
    ii = is_in_box(qs316,radecbox) # CS Find the sweep files we need to get all of points in the sweeps 
    if np.any(ii):
        sweeplist_in.append(fns)

print("Sweep List in:",sweeplist_in)
print("len(sweeplist_in):",len(sweeplist_in))

sweep_tables = []
for i in tqdm(range(len(sweeplist_in))):
    # CS Put sweep file objects in one table
    sweep = Table.read(sweeplist_in[i], memmap=True)
    mag_mask = (sweep["FLUX_R"] > rl19)
    sweep = sweep[mag_mask]
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
id1, id2, d2, d3 = c2.search_around_sky(c1, (1/3600)*u.deg) # CS 1" matching radius

qs316_match = qs316[id1] # CS quasar files matched objects
qs_crossmatch = sw_stack[id2] # CS Sweep files matched objects


print("qs_crossmatch len:",len(qs_crossmatch["RA"]))
print("qs316_match len:",len(qs316_match))
print("max(qs_crossmatch[RA]):",np.max(qs_crossmatch["RA"]))


# CS Everything looks good

# Mag cut
ii = (qs_crossmatch["FLUX_R"] > 19)
qs_cm_mc = qs_crossmatch[ii]
print(len(qs_cm_mc))
qs = qs_cm_mc # CS Shorter variable name

qs.write("qs.fits", overwrite=True)
"""

"""
qs = Table.read('qs.fits', memmap=True)
# CS So now qs_crossmatch is just quasars, but we have all
# CS of the info on them
# CS Now lets look for trends

g_z = []
g_r = []
z_w1 = []
r_w1 = []
r_w2 = []
r_w3 = []
r_w4 = []
flux_g = []
flux_r = []
flux_z = []
flux_i = []
flux_w1 = []
flux_w2 = []
flux_w3 = []
flux_w4 = []
pmra = []

fig, axs = plt.subplots(1, 1, figsize=(8,6))

for i, row in enumerate(qs):
    flux_g.append(row["FLUX_G"] / row["MW_TRANSMISSION_G"])
    flux_r.append(row["FLUX_R"] / row["MW_TRANSMISSION_R"])
    flux_z.append(row["FLUX_Z"] / row["MW_TRANSMISSION_Z"])
    flux_w1.append(row["FLUX_W1"] / row["MW_TRANSMISSION_W1"])
    flux_w2.append(row["FLUX_W2"] / row["MW_TRANSMISSION_W2"])
    flux_w3.append(row["FLUX_W3"] / row["MW_TRANSMISSION_W3"])
    flux_w4.append(row["FLUX_W4"] / row["MW_TRANSMISSION_W4"])
    pmra.append(row["PMRA"])

# CS Not really using parallax bc gaia parallax is crap at low values
# CS See X. Luri et al 2018 fig 8

flux_table = Table()

#plt.savefig("FluxTest.png")
for i in range(len(flux_g)):
    r_w1.append(flux_r[i]-flux_w1[i])
    z_w1.append(flux_z[i]-flux_w1[i])
    g_z.append(flux_g[i]-flux_z[i])
    g_r.append(flux_g[i]-flux_r[i])
    r_w2.append(flux_r[i]-flux_w2[i])
    r_w3.append(flux_r[i]-flux_w3[i])
    r_w4.append(flux_r[i]-flux_w4[i])


#plt.scatter(g_z,r_w1)
#plt.xlabel("g_z")
#plt.ylabel("r_w1")
#plt.savefig("gzVrw1")

gzm = np.mean(g_z)
gzs = np.std(g_z)

grOzw1 = []
for i in range(len(g_r)):
    grOzw1.append(g_r[i]/z_w1[i])


print("G_Z mean flux, G_Z stdev flux:",gzm,gzs)
print("G_R mean flux, G_R stdev flux:",np.mean(g_r),np.std(g_r))
print("R_W1 mean flux, R_W1 stdev flux:",np.mean(r_w1),np.std(r_w1))
print("Z_W1 mean flux, Z_W1 stdev flux:",np.mean(z_w1),np.std(z_w1))
print("R_W2 mean flux, R_W2 stdev flux:",np.mean(r_w2),np.std(r_w2))
print("Mean, stdev, G flux:", np.mean(flux_g),np.std(flux_g))
print("Mean, stdev, R flux:", np.mean(flux_r),np.std(flux_r))
print("Mean, stdev, Z flux:", np.mean(flux_z),np.std(flux_z))
print("Mean, stdev, w1 flux:", np.mean(flux_w1),np.std(flux_w1))
print("Mean, stdev, G_R/Z_W1:", np.mean(grOzw1),np.std(grOzw1))


print(len(qs))

#plt.scatter(qs["ANYMASK_G"][0],qs["ANYMASK_R"][0])
#plt.xlabel("ANYMASK_G")
#plt.ylabel("ANYMASK_G")
#plt.show()

# CS Lets look at saturation trends:

#CS recall bits are:
#000000 = 0
#000001 = 1
#000010 = 2
#000011 = 3

#2**2 = 4 bit 2
#2**3 = 8  bit 3
#2**4 = 16 bit 4

    
flag = 2**2 # CS Corresponds to saturation in g band in any exposure
ii = (qs["MASKBITS"] & flag) == 0
qsg = qs[ii]
    
flag = 2**3 # CS Corresponds to saturation in r band in any exposure
jj = (qs["MASKBITS"] & flag) == 0
qsr = qs[jj]


flag = 2**4 # CS Corresponds to saturation in z band
kk = (qs["MASKBITS"] & flag) == 0
qsz = qs[kk]

print("len(qsg):",len(qsg)) # CE Get out 190 so no saturated sources
print("len(qsr):",len(qsr))
print("len(qsz):",len(qsz))

# CS So none of these objects are saturated

flag = 2**5 # CS Corresponds to the source being a bright star
kk = (qs["FITBITS"] & flag) == 0
qsb = qs[kk] # CS These are sources that are NOT flagged as bright stars, which makes sense

flag = 2**6 # CS Corresponds to the source being a medium bright star
kk = (qs["FITBITS"] & flag) == 0
qsmb = qs[kk] # CS These are sources that are NOT flagged as medium bright stars, which makes sense
print("len(qsmb):",len(qsmb))


# CS So we only keep objects that aren't flagged as bright stars:
print("len(qsb):",len(qsb))
# CS Cool so none of them are flagged as bright stars, makes sense

flag = 2**7 # CS Corresponds to the source being a GAIA source
kk = (qs["FITBITS"] & flag) == 0
qsGS = qs[kk]

print("len(qsGS):",len(qsGS))
# CS Ok so only 2 of these are non GAIA sources
# CS So lets KEEP objects that are GAIA sources


print("Mean, stdv qs pmra:",np.mean(pmra),np.std(pmra))

"""

"""
FUNCTION SECTION

Things that signal quasar:
G_Z mean flux, G_Z stdev flux: -21.70186 43.12116
G_R mean flux, G_R stdev flux: -10.235811 17.815136
Z_W1 mean flux, Z_W1 stdev flux: -94.14673 255.5942
R_W1 mean flux, R_W1 stdev flux: -105.612785 275.46014, No change, not worth the comp time
R_W2 mean flux, R_W2 stdev flux: -189.35735 379.17133
Mean, stdev, G flux: 41.962494 27.61245
Mean, stdev, R flux: 52.198303 36.727577
Mean, stdev, Z flux: 63.66435 58.977596
Mean, stdev, w1 flux: 157.8111 303.04614
Mean, stdv qs para: 0.24795304 0.21279751
Mean, stdv qs pmra: 0.010143998 0.51280487
Mean, stdev, G_R/Z_W1: 0.047769446 5.487359


We take objects that are within 1 stdev of each flux measure

"""

def splendid_function(T_in):
    """ Applies a pre-determined color cut to
    to target quasars within a given astropy sweep

    Parameters
    ----------
    T_in : :class:`~Astropy Table`
        An astropy table that contains at minimum the same columns
        that are in a sweep file.

    Returns
    -------
    :class:`~numpy.ndarray`
        ``True`` for objects that are quasars, ``False`` for objects that are not quasars.
    """
    # CS Make an array to output
    array_out = []
    # CS Read in the table
    T = Table.read(str(T_in), memmap=True)
    # CS Matching Params:
    rw1m = -105.61 # CS Mean of r-w1
    rw1s = 275.46014 # CS Standard deviation of r-w1
    zw1m = -94.14673 # CS ect., ect.
    zw1s = 255.5942
    gzm = -21.70186 
    gzs = 43.12116
    grm = -10.235811
    grs = 17.815136
    rw2m = -189.35735
    rw2s = 379.17133
    mFg = 41.962494
    sFg = 27.61245
    mFr = 52.198303
    sFr = 36.727577
    mFz = 63.66435
    sFz = 58.977596
    mFw1 = 157.8111
    sFw1 = 303.04614
    grzw1m = 0.047769446 
    grzw1s = 5.487359
    
    para_cut = 2 # CS 2 mas = 500 pc
  
    mpmra = 0.010143998
    spmra = 0.51280487
    
    rl19 = float(10**(-(19-22.5)/2.5)) # CS Corresponds to flux =  r band mag of 19
    
    # CS Number of standard deviations to do our fit to
    nsd = 0.5
    
    # CS Calculate corrected fluxes in advance
    G = T["FLUX_G"] / T["MW_TRANSMISSION_G"]
    R = T["FLUX_R"] / T["MW_TRANSMISSION_R"]
    Z = T["FLUX_Z"] / T["MW_TRANSMISSION_Z"]
    W1 = T["FLUX_W1"] / T["MW_TRANSMISSION_W1"]
    
    master_mask = np.ones(len(T), dtype=bool)

    #print("Starting length:",len(T))
    #print("Keeping things that ARE flagged as GAIA sources")
    # CS Since most Quasars qre GAIA sources, lets keep just those
    master_mask &= (T["FITBITS"] & 2**7) != 0
        
    #print("Keeping sources that are NOT bright or medium bright stars")
    master_mask &= (T["FITBITS"] & 2**5) == 0 # CS Corresponds to the source being a bright star
    master_mask &= (T["FITBITS"] & 2**6) == 0 # CS Corresponds to the source being a medium bright star


    #print("Keeping sources that are NOT saturated in G, R, and Z")
    master_mask &= (T["MASKBITS"] & 2**2) == 0 # CS Corresponds to saturation in g band in any exposure
    master_mask &= (T["MASKBITS"] & 2**3) == 0 # CS Corresponds to saturation in r band in any exposure
    master_mask &= (T["MASKBITS"] & 2**4) == 0 # CS Corresponds to saturation in z band in any exposure

    #print("Now objects brigter with r < 19")
    master_mask &= (R > rl19) # CS Objects brighter than r 19

    #print("Now objects with G-Z within ", nsd," stdevs of G-Z mean")
    master_mask &= (np.abs((G-Z) - gzm) < nsd*gzs)
    
    #print("Now objects with R-W1 within ", nsd," stdevs of R-W1 mean")
    master_mask &= (np.abs((R-W1) - rw1m) < nsd*rw1s)
    
    #print("Now objects with G-R within ", nsd," stdevs of G-R mean")
    master_mask &= (np.abs((G-R) - grm) < nsd*grs)
    
    #print("Now objects with Z-W1 within ", nsd," stdevs of Z-W1 mean")
    master_mask &= (np.abs((Z-W1) - zw1m) < nsd*zw1s)
    
    #print("Now objects with (G-R)/(Z-W1)  within ", nsd," stdevs of (G-R)/(Z-W1) mean")
    master_mask &= (np.abs((G-R)/(Z-W1) - grzw1m) < nsd*grzw1s)

    
    #print("Now objects with fluxes within ", nsd," stdev of the qs mean")
    master_mask &= (np.abs(G-mFg) < nsd*sFg)
    master_mask &= (np.abs(R-mFr) < nsd*sFr)
    master_mask &= (np.abs(Z-mFz) < nsd*sFz)
    master_mask &= (np.abs(W1-mFw1) < nsd*sFw1)

    #print("Now objects with PMRAs within ", nsd," stdv of mean")
    master_mask &= (np.abs(T["PMRA"]-mpmra) < nsd*spmra)
    
    #print("Keeping objects that are further than 500pc")
    # CS This is around where I stop trusting them fully
    master_mask &= (T["PARALLAX"] < para_cut) # CS Keep negative parallaxes
    
    #print("Keeping objects that are TYPE = PSF")
    master_mask &= (T["TYPE"] == "PSF")
    
    T = T[master_mask]
    
    #print("splendid_function recovered:",len(T),"quasars")

    #print(master_mask[0:100])
    return master_mask
    
    
# CS args.Table = table input as python HW4.py mTable.fits
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("Table")
    parser.add_argument("--Description:", help = "Returns an array of True/False values, with True corresponding to quasars in the input astropy table that has, at minimum SDSS sweep file columns")
    args = parser.parse_args()
    splendid_function(args.Table)
    
    """
    CS Code testing block
    ii = splendid_function(args.Table)
    Ttest = Table.read('/d/scratch/ASTR5160/data/legacysurvey/dr9/south/sweep/9.0/sweep-170p025-180p030.fits', memmap=True)
    qout = Ttest[ii]
    qs316 = Table.read('qs.fits', memmap=True)
    
    qsra = np.array(qs316["RA"])*u.degree # CS 316 Quasars
    qsdec = np.array(qs316["DEC"])*u.degree
    c1 = SkyCoord(qsra,qsdec, frame='icrs')

    test_ra = np.array(Ttest["RA"])*u.degree # CS our combined sweep files
    test_dec = np.array(Ttest["DEC"])*u.degree # CS wtih coords < 3 deg of center
    c2 = SkyCoord(test_ra,test_dec,frame='icrs')
    id1, id2, d2, d3 = c2.search_around_sky(c1, (1/3600)*u.deg) # CS 1" matching radius

    qs316_match = qs316[id1] # CS quasar files matched objects
    print("Matched sources:",len(qs316_match))
    print("I got 68...")
    """

"""
CS Lets test

I get 190 r<19 quasars from these 4 sweep files
'/d/scratch/ASTR5160/data/legacysurvey/dr9/south/sweep/9.0/sweep-170p025-180p030.fits', My code gives: 68
'/d/scratch/ASTR5160/data/legacysurvey/dr9/south/sweep/9.0/sweep-170p030-180p035.fits', My code gives: 50
'/d/scratch/ASTR5160/data/legacysurvey/dr9/south/sweep/9.0/sweep-180p025-190p030.fits', My code gives: 59
'/d/scratch/ASTR5160/data/legacysurvey/dr9/south/sweep/9.0/sweep-180p030-190p035.fits', My code gives: 49
226, not bad
Let's see how many are "right"

1st sweep:
Now have: 68 Objects
Matched sources: 61
so ~90%!
"""