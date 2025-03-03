w# -*- coding: utf-8 -*-
"""
Spherical Caps

Chase L. Smith
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
        print("4-Vector representation is:", x,y,z,1)
        
#v4capr("05h00m00s")
#4-Vector representation is: -0.9659258262890682 0.258819045102521 0.0 1

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
        print("4-Vector representation for DEC:",dec," is:", f"{x:4.2f}",f"{y:4.2f}",f"{z:4.2f}",(1-math.sin(decc*math.pi/180)))
        
v4capd("36d00m00s")
#4-Vector representation is: 0.0 0.0 1.0 0.41221474770752686

if __name__ == "__main__":
    def v4caprd(ra,dec,radius):
        """
        This funciton to creates the vector-4 array
        that represents a circular field drawn
        on the surface of the sphere at a given ra and dec,
        that has some radius
        """
        c = SkyCoord(ra,dec, frame = 'icrs')
        c.representation_type = "cartesian"
        x = c.x
        y = c.y
        z = c.z
        print("4-Vector representation is:", x, y, z, 1-np.cos(radius*math.pi/180))

#v4caprd("05h00m00s","36d00m00s",1)
#4-Vector representation is: 0.20938900595583548 0.7814504087735196 0.5877852522924731 0.0001521504625372483

    
    
    
    
    