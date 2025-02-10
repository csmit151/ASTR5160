#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Week 4

Distances on the Sphere
"""

from astropy import units as u
from astropy.coordinates import SkyCoord
import numpy as np

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



#SepDot(c1.x,c1.y,c1.z,c2.x,c2.y,c2.z)
#Returns 0.880428183050399 rad

#Skycoord Seperation
#c1c2sep = c1.seperation(c2)
#print("C1 and C2 sep:",c1c2sep.radian)