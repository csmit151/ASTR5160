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
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
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

sweep1 = "/d/scratch/ASTR5160/data/legacysurvey/dr9/south/sweep/9.0/sweep-170p025-180p030.fits"
sweep2 = "/d/scratch/ASTR5160/data/legacysurvey/dr9/south/sweep/9.0/sweep-170p030-180p035.fits"
sweep3 = "/d/scratch/ASTR5160/data/legacysurvey/dr9/south/sweep/9.0/sweep-180p025-190p030.fits"
sweep4 = "/d/scratch/ASTR5160/data/legacysurvey/dr9/south/sweep/9.0/sweep-180p030-190p035.fits"


c1 = SkyCoord(stars["RA"]*u.deg,stars["DEC"]*u.deg,frame='icrs')
c2 = SkyCoord(qs["RA"]*u.deg,qs["DEC"]*u.deg,frame='icrs')
id1,id2,d2,d3 = c1.search_around_sky(c2,(0.5/3.6)*u.deg)


#plt.scatter(stars["RA"][id1],stars["DEC"][id1])
#plt.scatter(qs["RA"][id2],qs["DEC"][id2],color="r")

gall = []
zall = []
rall = []
w1all = []
swra = []
swdec = []
sw1ra = []
sw1dec = []
sw2ra = []
sw2dec = []
sw3ra = []
sw3dec = []
sw4ra = []
sw4dec = []

with fits.open(sweep1) as hdul:
    data = hdul[1].data
    sw1ra.append(data["RA"])
    sw1dec.append(data["DEC"])
    g1 = data["FLUX_G"]/data["MW_TRANSMISSION_G"]
    r1 = data["FLUX_R"]/data["MW_TRANSMISSION_R"]
    z1 = data["FLUX_Z"]/data["MW_TRANSMISSION_Z"]
    w1_1 = data["FLUX_W1"]/data["MW_TRANSMISSION_W1"]
    w2_1 = data["FLUX_W2"]/data["MW_TRANSMISSION_W2"]
    
with fits.open(sweep2) as hdul:
    data = hdul[1].data
    sw2ra.append(data["RA"])
    sw2dec.append(data["DEC"])
    g2 = data["FLUX_G"]/data["MW_TRANSMISSION_G"]
    r2 = data["FLUX_R"]/data["MW_TRANSMISSION_R"]
    z2 = data["FLUX_Z"]/data["MW_TRANSMISSION_Z"]
    w1_2 = data["FLUX_W1"]/data["MW_TRANSMISSION_W1"]
    w2_2 = data["FLUX_W2"]/data["MW_TRANSMISSION_W2"]

with fits.open(sweep3) as hdul:
    data = hdul[1].data
    sw3ra.append(data["RA"])
    sw3dec.append(data["DEC"])
    g3 = data["FLUX_G"]/data["MW_TRANSMISSION_G"]
    r3 = data["FLUX_R"]/data["MW_TRANSMISSION_R"]
    z3 = data["FLUX_Z"]/data["MW_TRANSMISSION_Z"]
    w1_3 = data["FLUX_W1"]/data["MW_TRANSMISSION_W1"]
    w2_3 = data["FLUX_W2"]/data["MW_TRANSMISSION_W2"]

with fits.open(sweep4) as hdul:
    data = hdul[1].data
    sw4ra.append(data["RA"])
    sw4dec.append(data["DEC"])
    g4 = data["FLUX_G"]/data["MW_TRANSMISSION_G"]
    r4 = data["FLUX_R"]/data["MW_TRANSMISSION_R"]
    z4 = data["FLUX_Z"]/data["MW_TRANSMISSION_Z"]
    w1_4 = data["FLUX_W1"]/data["MW_TRANSMISSION_W1"]
    w2_4 = data["FLUX_W2"]/data["MW_TRANSMISSION_W2"]

sw1raz = sw1ra[0]
sw2raz = sw2ra[0]
sw3raz = sw3ra[0]
sw4raz = sw4ra[0]
for i in range(len(sw1raz)):
    swra.append(float(sw1raz[i]))
for i in range(len(sw2raz)):
    swra.append(float(sw2raz[i]))
for i in range(len(sw3raz)):
    swra.append(float(sw3raz[i]))
for i in range(len(sw4raz)):
    swra.append(float(sw4raz[i]))

sw1decz = sw1dec[0]
sw2decz = sw2dec[0]
sw3decz = sw3dec[0]
sw4decz = sw4dec[0]

for i in range(len(sw1decz)):
    swdec.append(float(sw1decz[i]))
for i in range(len(sw2decz)):
    swdec.append(float(sw2decz[i]))
for i in range(len(sw3decz)):
    swdec.append(float(sw3decz[i]))
for i in range(len(sw4decz)):
    swdec.append(float(sw4decz[i]))


for i in range(len(g1)):
    gall.append(float(g1[i]))
for i in range(len(g2)):
    gall.append(float(g2[i]))
for i in range(len(g3)):
    gall.append(float(g3[i]))
for i in range(len(g4)):
    gall.append(float(g4[i]))

for i in range(len(z1)):
    zall.append(float(z1[i]))
for i in range(len(z2)):
    zall.append(float(z2[i]))
for i in range(len(z3)):
    zall.append(float(z3[i]))
for i in range(len(z4)):
    zall.append(float(z4[i]))

for i in range(len(r1)):
    rall.append(float(r1[i]))
for i in range(len(r2)):
    rall.append(float(r2[i]))
for i in range(len(r3)):
    rall.append(float(r3[i]))
for i in range(len(r4)):
    rall.append(float(r4[i]))

for i in range(len(w1_1)):
    w1all.append(float(w1_1[i]))
for i in range(len(w1_2)):
    w1all.append(float(w1_2[i]))
for i in range(len(w1_3)):
    w1all.append(float(w1_3[i]))
for i in range(len(w1_4)):
    w1all.append(float(w1_4[i]))

print(len(gall))
print(len(swra),len(swdec))

swracor = []
swdeccor = []
for i in range(len(swra)):
    if (swra[i] > 177) & (swra[i] < 183) & (swdec[i] > 27) & (swdec[i] < 33):
        swracor.append(swra[i])
        swdeccor.append(swdec[i])
            
print(len(stars["RA"][id1]))

gfinal = []
zfinal = []
rfinal = []
w1final = []

# CS now convert the fluxes for each confimed qs and star to mags

print(len(id1))
print(len(swracor))


c3 = SkyCoord(stars["RA"][id1]*u.deg,stars["DEC"][id1]*u.deg,frame='icrs')
c4 = SkyCoord(swracor*u.deg,swdeccor*u.deg,frame='icrs')
id3,id4,d3,d4 = c3.search_around_sky(c4,(0.5)*u.deg)


for i in tqdm(range(len(id3))):
    gfinal.append(22.5-2.5*np.log10(float(gall[id3[i]])))        
    zfinal.append(22.5-2.5*np.log10(float(zall[id3[i]])))
    rfinal.append(22.5-2.5*np.log10(float(rall[id3[i]])))
    w1final.append(22.5-2.5*np.log10(float(w1all[id3[i]])))

# CS plot r-W1 vs g-z

r_w1 = []
g_z = []
for i in range(len(rfinal)):
    r_w1.append(rfinal[i]-w1final[i])
    g_z.append(gfinal[i]-zfinal[i])
    
plt.scatter(g_z,r_w1)
plt.xlabel("G-Z")
plt.ylabel("R-W1")
plt.show()


#if __name__ == "__main__":
    #print("Coord matched star sweep files:")
    #whichSweepSouth(stars)
    #print("Coord matched qs sweep files:")
    #whichSweepSouth(qs)
    
    
