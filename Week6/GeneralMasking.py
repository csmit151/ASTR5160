#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
General Masking

Chase L. Smith
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

#CS Note that if theta > theta_c (as in 1-h=1-cos(theta_c))
#then the point lise outside of the cap, and visa versa for inside the cap



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


if __name__ == "__main__":
    def v4capr(ra):
        """
        This funciton to creates the vector-4 array
        x,y,z,1-cos(theta) for the spherical cap bounded by some input RA
        """
        
        c = SkyCoord(ra,'00d00m00s', frame = 'icrs')
        ra = c.ra
        dec = c.dec
        cc = SkyCoord(ra+(90*u.deg),dec, frame = 'icrs')
        cc.representation_type = "cartesian"
        x = cc.x
        y = cc.y
        z = cc.z
        return(x,y,z,1)
        

if __name__ == "__main__":
    def v4capd(dec):
        """
        This funciton to creates the vector-4 array
        for the spherical cap bounded by some input DEC
        """
        
        c = SkyCoord('00h00m00s',dec, frame = 'icrs')
        decc = c.dec.degree
        c.representation_type = "cartesian"
        c2 = SkyCoord('00h00m00s','90d00m00s', frame = 'icrs')
        c2.representation_type = "cartesian"
        x = c2.x
        y = c2.y
        z = c2.z
        return(x,y,z,(1-math.sin(decc*math.pi/180))) 
        
#CS Create a polygon in a Mangle file consisting of 4 caps
#RA 5h, 6h, dec 30 deg, 40 deg

x1, y1, z1, h1 = v4capr("05h00m00s")
x2, y2, z2, h2 = v4capr("06h00m00s")
x3, y3, z3, h3 = v4capd("30d00m00s")
x4, y4, z4, h4 = v4capd("40d00m00s")

mf = open("mask.ply", "w")
mf.write("1 polygons"+"\n")
mf.write("polygon 1 ( 4 caps, 1 weight, 0 pixel, 0 str):\n")
mf.write("  "+str(x1)+" "+str(y1)+" "+str(z1)+" "+str(h1)+"\n")
mf.write("  "+str(x2)+" "+str(y2)+" "+str(z2)+" "+str(h2)+"\n")
mf.write("  "+str(x3)+" "+str(y3)+" "+str(z3)+" "+str(h2)+"\n")
mf.write("  "+str(x4)+" "+str(y4)+" "+str(z4)+" "+str(h4)+"\n")
mf.close()

m = pymangle.Mangle("mask.ply")
ra = np.array([47.3, 152.7, 23.3, 280.4])
dec = np.array([-11.2, 12.2, 88.7, -39.2])
good = m.contains(ra, dec)

ra_rand, dec_rand = m.genrand(10000)

#CS Area of lat-lon rec:
def Alatlon(ra1,ra2,dec1,dec2):
    """
    CS This function returns the area of
    a lat-lon rectangle in str given input ra in radians
    and dec in radians
    """
    #print("Area:",(ra2-ra1)*(math.sin(dec2)-math.sin(dec1)))
    return (ra2-ra1)*(math.sin(dec2)-math.sin(dec1))
    

#Calculate correct str area for my polygon
Alatlon((5*15*np.pi/180),(6*15*np.pi/180),(30*np.pi/180),(40*np.pi/180))
al1 = Alatlon((5*15*np.pi/180),(6*15*np.pi/180),(30*np.pi/180),(40*np.pi/180))

#CS Area: 0.03738170880123988

#CS Add a second polygon to the file for a field bounded in RA
#by 11h and 12h and dec in 60d 70d
x5, y5, z5, h5 = v4capr("11h00m00s")
x6, y6, z6, h6 = v4capr("12h00m00s")
x7, y7, z7, h7 = v4capd("60d00m00s")
x8, y8, z8, h8 = v4capd("70d00m00s")
al2 = Alatlon((11*15*np.pi/180),(12*15*np.pi/180),(60*np.pi/180),(70*np.pi/180))


mf = open("mask.ply", "w")
mf.write("2 polygons"+"\n")
mf.write("polygon 1 ( 4 caps, 1 weight, 0 pixel, "+str(al1)+" str):\n")
mf.write("  "+str(x1)+" "+str(y1)+" "+str(z1)+" "+str(h1)+"\n")
mf.write("  "+str(x2)+" "+str(y2)+" "+str(z2)+" "+str(h2)+"\n")
mf.write("  "+str(x3)+" "+str(y3)+" "+str(z3)+" "+str(h2)+"\n")
mf.write("  "+str(x4)+" "+str(y4)+" "+str(z4)+" "+str(h4)+"\n")
mf.write("polygon 2 ( 4 caps, 1 weight, 0 pixel, "+str(al2)+" str):\n")
mf.write("  "+str(x5)+" "+str(y5)+" "+str(z5)+" "+str(h5)+"\n")
mf.write("  "+str(x6)+" "+str(y6)+" "+str(z6)+" "+str(h6)+"\n")
mf.write("  "+str(x7)+" "+str(y7)+" "+str(z7)+" "+str(h7)+"\n")
mf.write("  "+str(x8)+" "+str(y8)+" "+str(z8)+" "+str(h8)+"\n")

mf.close()

