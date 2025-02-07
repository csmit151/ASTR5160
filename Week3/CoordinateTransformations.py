# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 16:35:50 2025

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

#Chase L. Smith
#Coordinate Transformation

"""
Convert some RA and DEC to Cartesian coordinates
"""
c = SkyCoord('02h37m23s','-37d17m32s', frame = 'icrs')
print("Converted Cooridantes are:",c,"Deg")

c.representation_type = "cartesian"
print("SkyCoord c in cart:",c)

#Check to see if this agrees with equations
RAdeg = 39.34583333
DECdeg = -37.29222222

x = math.cos(RAdeg*(3.14)/180)*math.cos(DECdeg*(3.14)/180)
y = math.sin(RAdeg*(3.14)/180)*math.cos(DECdeg*(3.14)/180)
z = math.sin(DECdeg)

print("From Adam's eq:",x,y,z)
#So everything checks out

#calculate the RA and Dec of the center of the galaxy
cgal = SkyCoord(frame='galactic', l=0*u.deg, b=0*u.deg)
cgalradec = cgal.transform_to("icrs")

print("Center of Galaxy:",cgalradec)

"""
Cent gal is 266.40, -28.96
Or: 17:45, -28.96

The center of the galaxy is in the constelation
Sagitatius(SagA* is 17h45m40s, -29 so this checks out),
and it is located roughly near the edge of the constelation
"""

#For laramie, DEC=40N, plot how l,b change through the year at zenith
#(i.e. directly above your head)
for i in range(0,360,1):
    clar = SkyCoord(float(i)*u.deg,40*u.deg,frame="icrs") #0 deg is the vernal equinox
    clargalactic = clar.transform_to("galactic")
    #print("Laramie zenith at",i,"days from the Vernal Equinox:",clargalactic)
    plt.scatter(clargalactic.l,clargalactic.b, color='blue')
    plt.xlabel("l Galactic Coordinate")
    plt.ylabel("b Galactic Coordinate")


