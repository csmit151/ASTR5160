#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Week 4

#Map Projections
Chase L. Smith
"""


#generate a random set of 10,000 points on a sphere
import matplotlib.pyplot as plt
import numpy as np
from numpy.random import random

N=10**4

def PonSph(N):
    """
    This function generates N
    random points on the surface of a sphere
    with coordinates ra, dec in radians
    
    This funciton populates the sphere equally in area
    on an XY grid
    """
    
    #Generate ra and dec
    ra=2*np.pi*(random(N)-0.5)
    dec = np.arcsin(1.-random(N)*2)
    #Plot the points
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(ra,dec, marker='o',color='b',s=0.7,alpha=0.5)
    ax.set_xlabel("RA in Radians")
    ax.set_ylabel("DEC in Radians")
    ax.grid(color='k',linestyle='solid',linewidth=0.5)
    #Show Fig
    fig.show()

#PonSph(N)
"""
There are more points near the equator of the sphere than the poles
"""


#Now plot points in an Aitoff projection

def PonSphA(N):
    """
    This function generates N
    random points on the surface of a sphere
    with coordinates ra, dec in radians in an Aitoff
    projection
    
    This funciton also populates the sphere equally in area
    
    This function includes a thick, blue, dashed axis grid
    """
    
    #Generate ra and dec
    ra=2*np.pi*(random(N)-0.5)
    dec = np.arcsin(1.-random(N)*2)
    #Plot the points
    fig = plt.figure()
    ax = fig.add_subplot(111,projection="aitoff")
    ax.scatter(ra,dec, marker='o',color='orange',s=0.7,alpha=0.5)
    xlab = ['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h']
    ax.set_xticklabels(xlab, weight=800)
    ax.grid(color='b',linestyle='dashed',linewidth=2)
    #Show Fig
    #fig.show()
    #Return an array of ra and dec points
    coordsra = []
    coordsdec = []
    for i in range(len(ra)):
        coordsra.append(ra[i])
        coordsdec.append(dec[i])
    return coordsra, coordsdec

#PonSphA(N)

#Now with a Lambert projection
def PonSphL(N):
    """
    This function generates N
    random points on the surface of a sphere
    with coordinates ra, dec in radians in an Lambert
    projection
    
    This funciton also populates the sphere equally in area
    
    This function includes a thin, blue, solid axis grid
    """
    
    #Generate ra and dec
    ra=2*np.pi*(random(N)-0.5)
    dec = np.arcsin(1.-random(N)*2)
    #Plot the points
    fig = plt.figure()
    ax = fig.add_subplot(111,projection="lambert")
    ax.scatter(ra,dec, marker='o',color='green',s=0.7,alpha=0.5)
    xlab = ['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h']
    ax.set_xticklabels(xlab, weight=800)
    ax.grid(color='b',linestyle='solid',linewidth=1)
    #Show Fig
    fig.show()

PonSphL(N)



