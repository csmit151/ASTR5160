#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
HW2

@author: csmit151
"""

# CS Import relevant packages
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.coordinates import AltAz
from astropy.time import Time
from astropy import units as u
import math
import matplotlib.patches as mpatches
from Week4.MapProjections import PonSphA
import argparse


# CS Use argparse to take an input file directory to save images
# CS from the command line
parser = argparse.ArgumentParser()
parser.add_argument("dname") # CS destination file path name
args = parser.parse_args()
dpath = args.dname

# CS Writes image to a location from cmd line in the form:
# CS python python HW2.py /d/ori1/csmit/AdamClass/imtest/

# CS Begin first function

def AreaLatLon(ramin,ramax,decmin,decmax):
    """Finds the area of a lat-lon rectangular feild in the sky that has
    corners in RA and DEC, and plots the rectangle in Aitoff projection
    and labels it with the rectangles area

    Parameters
    ----------
    ramin : :class:`int`
        The minimum RA value of the box (degrees).
    ramax : :class:`int`
        The maxmimum  RA value of the box (degrees).
    decmin : :class:`int`
        The minimum  DEC value of the box (degrees).
    decmax : :class:`int`
        The maxmimum  DEC value of the box (degrees).

    Returns
    -------
    :class:`float`
        Area of the lat-lon rectangle in square degrees

    """
    # CS Calculate the area of a rectangular feild w/ class 9 notes
    Area = (180/math.pi)*(ramax-ramin)*(math.sin(decmax*math.pi/180)-math.sin(decmin*math.pi/180))
    print("Area of lat-lon rectangle in square deg:", Area)
    
    # CS Make a list of ra and dec points that are along the sides of our rectangle
    s1 = np.arange(ramin,ramax,0.01) # CS points along our "sides"
    s2 = np.arange(decmin,decmax,0.01)
    dminside = np.full(len(s1),decmin) # CS each ra "side" needs to be at both decmin and decmas, and visa versa for dec sides
    dmaxside = np.full(len(s1),decmax)
    raminside = np.full(len(s2),ramin)
    ramaxside = np.full(len(s2),ramax)
    
    # CS Choose a random color for the lines for easier viewing
    colors = ['b','g','r','c','m','y','orange']
    cnum = np.random.randint(0,len(colors))
    # CS plot four sides, each anchored on the opposing coords min/max value
    rf = math.pi/180 # CS factor for converting deg to radians
    # CS Begin plotting
    fig = plt.figure()
    # CS Make sure we're in Aitoff projection
    ax = fig.add_subplot(111,projection="aitoff")
    # CS Plot the four sides of our lat-lon rectangle, making sure we're in radians
    ax.scatter(s1*rf,dminside*rf, marker='o',color=colors[cnum],s=0.7,alpha=0.5)
    ax.scatter(s1*rf,dmaxside*rf, marker='o',color=colors[cnum],s=0.7,alpha=0.5)
    ax.scatter(raminside*rf,s2*rf, marker='o',color=colors[cnum],s=0.7,alpha=0.5)
    ax.scatter(ramaxside*rf,s2*rf, marker='o',color=colors[cnum],s=0.7,alpha=0.5)
    # CS Now label the region with its area. 
    Alabel = mpatches.Patch(facecolor=colors[cnum], label = "Area: "+str(np.round(Area, decimals=2))+" square deg", linewidth = 0.5, edgecolor = 'black')
    legend = plt.legend(handles=[Alabel], loc = 9, bbox_to_anchor = (0.95,1), fontsize = 12, fancybox = False, edgecolor = "none", facecolor = "none")
    
    # CS Set up hour labels every 2hrs
    xlab = ['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h']
    ax.set_xticklabels(xlab, weight=800)
    ax.grid(color='b',linestyle='dashed',linewidth=2)
    # CS show Fig
    fig.show()
    plt.savefig(str(dpath)+"ALL_"+str(ramin)+str(ramax)+str(decmin)+str(decmax)+".png")
    # CS return the "Area" value
    return Area

    
"""
CS For #1: Area of lat-lon rectangle in square deg: 20626.480624709635
Given 0,360,0,90 so this makes sense given that that is exactly 1/2 of 
the total area of a sphere in square degrees, 41252.96
"""

# CS Begin second function

def PopLatLon(ramin,ramax,decmin,decmax):
    """Randomly populates a lat-lon rectangular feild drawn on the surface
    of a sphere with 10**4 points. The feild is populated randomly in area. 

    Parameters
    ----------
    ramin : :class:`int`
        The minimum RA value of the box (degrees).
    ramax : :class:`int`
        The maxmimum  RA value of the box (degrees).
    decmin : :class:`int`
        The minimum  DEC value of the box (degrees).
    decmax : :class:`int`
        The maxmimum  DEC value of the box (degrees).

    Returns
    -------
    :class:`array`
        A set of (ra,dec) coordinates that lie within
        the lat-lon rectangular feild. 
    """
    
    # CS first start by populating the entire sphere randomly in area
    # CS Import PonSphA from Week4
    # CS This function generates Nrandom points on the surface of a sphere
    # CS with coordinates ra, dec in radians in an Aitoff projection
    # CS This funciton also populates the sphere equally in area
    coordsra, coordsdec = PonSphA(10**4)
    # CS Now we need to just return the set of ra and dec within our rectangle
    rf = math.pi/180 # CS factor for converting deg to radians
    ii = []
    
    # CS Record index of which points are in our rectangle
    for k in range(len(coordsra)):
        if (coordsra[k] > ramin*rf) and (coordsra[k] < ramax*rf) and (coordsdec[k] > decmin*rf) and (coordsdec[k] < decmax*rf):
            ii.append(k)
    
    AreaLatLon(ramin,ramax,decmin,decmax)
    # CS Begin plotting
    fig = plt.figure()
    # CS Make sure we're in Aitoff projection
    ax = fig.add_subplot(111,projection="aitoff")
    # CS Plot the points within our rectangle
    for i in range(len(ii)):
        ax.scatter(coordsra[ii[i]],coordsdec[ii[i]], marker='o',color='r',s=5)
    #CS Show fig
    fig.show()
    plt.savefig(str(dpath)+"PLL_"+str(ramin)+str(ramax)+str(decmin)+str(decmax)+".png")
    
    #print("% filled:",(float(len(ii))/10**4)*100)
    
if __name__ == "__main__":
    AreaLatLon(15,60,10,40)
    AreaLatLon(15,60,30,60)
    AreaLatLon(15,60,50,80)
    PopLatLon(15,60,10,40)
    #print(1209.586/41252.961*100, "%")
    
    
"""
CS For #2: for:
    AreaLatLon(15,60,10,40)
    PopLatLon(15,60,10,40)
I get an area from AreaLatLon of Area of lat-lon rectangle in square deg:', 1209.5869256049684
The total area of a sphere is 41252.961 square deg
So this gives us a rectangle covering 2.9321192241206635 % of the sphere
Comparing with the % of random points that land within our rectangle, we find:
'% filled:', 2.91 (Note that this number fluctuates slightly, a few runs are below)
2.86
3.05
3.01

CS So we confirm that the lat-lon rectangle that we are populating
contains roughly the correct number of random points

"""
    
    