#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 11:50:28 2025

@author: csmit151
"""


#CS Import relevant packages
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.coordinates import AltAz
from astropy.time import Time
from astropy import units as u
import math
import pymangle


if __name__ == "__main__":
    def v4caprd(ra,dec,radius):
        """
        CS This funciton to creates the vector-4 array
        that represents a circular field drawn
        on the surface of the sphere at a given ra and dec,
        that has some radius
        """
        #CS make a SkyCoord from input ra and dec
        c = SkyCoord(ra,dec, frame = 'icrs')
        #CS extract cartesian coordinates
        c.representation_type = "cartesian"
        x = c.x
        y = c.y
        z = c.z
        #print(x, y, z, 1-np.cos(radius*math.pi/180))
        return x, y, z, 1-np.cos(radius*math.pi/180)

x1, y1, z1, h1 = v4caprd("76d00m00s","36d00m00s",5)
x2, y2, z2, h2 = v4caprd("75d00m00s","35d00m00s",5)


"""
CS This returns:
4-Vector representation is: 0.19571892485153297 0.7849857321266633 0.5877852522924731 0.003805301908254455
4-Vector representation is: 0.21201214989665462 0.7912401152362238 0.573576436351046 0.003805301908254455
"""

#CS Create a file contating these two caps as intersection.ply 
#in the mangle format


mf = open("intersection.ply", "w")
mf.write("1 polygons"+"\n")
mf.write("polygon 1 ( 2 caps, 1 weight, 0 pixel, 0 str):\n")
mf.write("  "+str(x1)+" "+str(y1)+" "+str(z1)+" "+str(h1))
mf.write("\n")
mf.write("  "+str(x2)+" "+str(y2)+" "+str(z2)+" "+str(h2))
mf.close()

#CS Create a second Mangle file contaning these two caps as diffrent polygons
mf = open("bothcaps.ply", "w")
mf.write("2 polygons"+"\n")
mf.write("polygon 1 ( 1 caps, 1 weight, 0 pixel, 0 str):\n")
mf.write("  "+str(x1)+" "+str(y1)+" "+str(z1)+" "+str(h1))
mf.write("\n")
mf.write("polygon 2 ( 1 caps, 1 weight, 0 pixel, 0 str):\n")
mf.write("  "+str(x2)+" "+str(y2)+" "+str(z2)+" "+str(h2))
mf.close()


#CS import and read in each mask
minter=pymangle.Mangle("intersection.ply")
mboth = pymangle.Mangle("bothcaps.ply","w")

#CS use genrand to fill each mask with 10,000 rand points
ra_rand_in, dec_rand_in = minter.genrand(10000)
ra_rand_b, dec_rand_b = mboth.genrand(10000)

#CS Now plot the random points corresponding to each mask
#plt.scatter(ra_rand_in,dec_rand_in,color="red", s=10)
#plt.scatter(ra_rand_b,dec_rand_b,color="blue", s=10)
#plt.savefig("/d/ori1/csmit/AdamClass/Mangle.png",dpi=300)
#plt.show()


#CS Flip the sign of the constraint on cap 1
mf = open("intersectionflip.ply", "w")
mf.write("1 polygons"+"\n")
mf.write("polygon 1 ( 2 caps, 1 weight, 0 pixel, 0 str):\n")
mf.write("  "+str(x1)+" "+str(y1)+" "+str(z1)+" -"+str(h1))
mf.write("\n")
mf.write("  "+str(x2)+" "+str(y2)+" "+str(z2)+" "+str(h2))
mf.close()

mflip1 = pymangle.Mangle("intersectionflip.ply")
ra_rand_inf, dec_rand_inf = mflip1.genrand(10000)

#CS plot minter and mflip1 on one plot in diffrent colors
plt.scatter(ra_rand_in,dec_rand_in,color="red", s=10)
plt.scatter(ra_rand_inf,dec_rand_inf,color="green", s=10)
plt.savefig("/d/ori1/csmit/AdamClass/MangleFlip.png",dpi=300)
plt.show()
