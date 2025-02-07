# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 13:20:06 2025

@author: UW-User
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy import units as u
import math
from astropy.time import Time
from astropy.coordinates import EarthLocation
import matplotlib.patches as mpatches


#Dust maps

import dustmaps
from dustmaps.config import config
from dustmaps.sfd import SFDQuery

ra, dec = '12h42m30s', '+41d12m00s' 
c = SkyCoord(ra, dec).galactic 

#Obtain the reddening at this ra dec position from the dust maps
dustdir = "/d/scratch/ASTR5160/data/dust/v0_1/maps"
config["data_dir"] = dustdir
sfd = SFDQuery()

print("sfd(c):",sfd(c))    #0.016697595

#Now we can obtain the reddening without first converting to Galactic coordinates

ebv = sfd(c)
ugriz = np.array([4.239 ,3.303 , 2.285 , 1.698 , 1.263])
A = ebv*ugriz

#with:
ra1, dec1 = "246.933", "40.795" #and 
ra2, dec2 = "236.562", "2.440"  
#are both quasars near a redshift of
z = 1.08

#From SDSS
#Quasar 1:
u1 = 18.82
g1 = 18.81
r1 = 18.73
i1 = 18.82
z1 = 18.90

#Quasar 2:
u2 = 19.37
g2 = 19.10
r2 = 18.79
i2 = 18.73
z2 = 18.63

#Now plot g-r and r-i
plt.scatter(g1-r1,r1-i1, marker="*", color='r', s = 80)
plt.scatter(g2-r2,r2-i2, marker="*", color='b', s = 80)

"""
These quasars don't have the same color,
thought they should as they are both at the same redshift.
The diffrence is dust reddening

Correcting the Quasar's color for dust extinction and replotting:
"""


c1 = SkyCoord(ra1,dec1).galactic
ebv1 = sfd(c1)
ugriz1 = np.array([u1, g1, r1, i1, z1])
A1 = ebv1*ugriz1

c2 = SkyCoord(ra2,dec2).galactic
ebv2 = sfd(c2)
ugriz2 = np.array([u2, g2, r2, i2, z2])
A2 = ebv2*ugriz2


#dust corrected ugriz values
ugriz1cor = np.add(ugriz1,A1)
ugriz2cor = np.add(ugriz2,A2)

#Now replot
plt.xlabel("g-r", fontsize=14)
plt.ylabel("r-i", fontsize=14)
plt.xlim(0,0.4)
plt.ylim(-0.3,0.3)
Q1 = mpatches.Patch(facecolor='red', label='Q1', linewidth = 0.5, edgecolor = 'black')
Q2 = mpatches.Patch(facecolor='blue', label='Q2', linewidth = 0.5, edgecolor = 'black')
Q1c = mpatches.Patch(facecolor='darkred', label='Q1 Dust Cor.', linewidth = 0.5, edgecolor = 'black')
Q2c = mpatches.Patch(facecolor='darkblue', label='Q2 Dust Cor.', linewidth = 0.5, edgecolor = 'black')
legend = plt.legend(handles=[Q1,Q2,Q1c,Q2c], loc = 0, bbox_to_anchor = (1,1), fontsize = 12, fancybox = False, edgecolor = "none", facecolor = "none") #loc code 9 is upper center, 2 upper left



#Visualize the dust in the region of each quasar

#2d mesh grid 100x100 array
#Centered on 236.6, 2.4
cenra = (236.6-5) #We subtract 5 in order to center our array on 236.6
cendec = (2.4-5)

meshrange = np.arange(0,10,0.1)
ragrid = []
decgrid = []

for i in range(len(meshrange)):
    ragrid.append(float(cenra)+float(meshrange[i]))
    decgrid.append(float(cendec)+float(meshrange[i]))

mg1 = np.meshgrid(ragrid,decgrid)


#Centered on 246.9, 40.8
cenra2 = (246.9-6.5) #We subtract 6.5
cendec2 = (40.8-5) #and 5 like before for dec

meshrangedec2 = np.arange(0,10,0.1) #0.1 bins for dec
meshrangera2 = np.arange(0,13,0.13) #0.13 bins for RA

ragrid2 = []
decgrid2 = []

for i in range(len(meshrangera2)):
    ragrid2.append(float(cenra2)+float(meshrangera2[i]))
    
for i in range(len(meshrangedec2)):
    decgrid2.append(float(cendec2)+float(meshrangedec2[i]))

mg2 = np.meshgrid(ragrid2,decgrid2)


