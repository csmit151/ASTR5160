#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Week 4

Distances on the Sphere
"""

from astropy import units as u
from astropy.coordinates import SkyCoord
import numpy as np
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

#Use astropy to convert (263.75, -17.9) and (20h24m59.9s,10d6m0s) to cartesian

c1 = SkyCoord(263.75*u.deg,-17.9*u.deg,frame='icrs')
c2 = SkyCoord("20h24m59.9s","10d6m0s",frame='icrs')



c1.representation_type = "cartesian"
c2.representation_type = "cartesian"
print("C1 cartesian:",c1)
print("C2 cartesian:",c2)

#Dot product:
def SepDot(x1,y1,z1,x2,y2,z2):
    """
    This function takes 6 input coordinates for two points on the sky
    in cartesian coordinates and retuns the Z_ang between them
    """
    print(np.arccos((x1*x2+y1*y2+z1*z2)/((x1**2+y1**2+z1**2)**(1/2)*(x2**2+y2**2+z2**2)**(1/2))))



SepDot(c1.x,c1.y,c1.z,c2.x,c2.y,c2.z)
#Returns 0.880428183050399 rad

#Skycoord Seperation
sep = c1.separation(c2)
print("C1 and C2 sep:",sep.radian)
#C1 and C2 sep: 0.8804281830503986

#use np.random.random to populate the area of the sky between ra=2-3hr and dec -2 to 2deg
#with two diffrent sets of 100 random points

#using arrays instead of loops:
a1 = [2]*100+np.array(np.random.random(100))
d1 = [-2]*100+4*np.array(np.random.random(100))

a2 = [2]*100+np.array(np.random.random(100))
d2 = [-2]*100+4*np.array(np.random.random(100))

#Plot both sets with diffrent colors/symbols
plt.scatter(a1,d1,color='red', marker="*")
plt.scatter(a2,d2,color='blue', marker="d")
plt.xlabel("RA (Hours)", fontsize=16)
plt.ylabel("DEC (deg)", fontsize=16)
rd1 = mpatches.Patch(facecolor='red', label='RA/DEC 1', linewidth = 0.5, edgecolor = 'black')
rd2 = mpatches.Patch(facecolor='blue', label='RA/DEC 2', linewidth = 0.5, edgecolor = 'black')
legend = plt.legend(handles=[rd1,rd2], loc = 0, bbox_to_anchor = (1,1), fontsize = 12, fancybox = False, edgecolor = "none", facecolor = "none") #loc code 9 is upper center, 2 upper left

#which points are within 10' of eachother?
crd1 = SkyCoord(a1*u.hr,d1*u.deg,frame='icrs')
crd2 = SkyCoord(a2*u.hr,d2*u.deg,frame='icrs')

id1, id2, dis1, dis2 = crd1.search_around_sky(crd2,1/6*u.deg)

#print(id1,id2,d1,d2)
#overplot there points in a new color
plt.scatter(a1[id1],d1[id1],color='green', marker="o",alpha=0.5, s=50)

