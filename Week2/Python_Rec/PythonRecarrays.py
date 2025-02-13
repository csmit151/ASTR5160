# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 15:51:03 2025

@author: UW-User
"""
#Python Recarrays

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table

#Python Recarrays Tasks
from astropy.table import Table

#Import stuc.fits file
objs = Table.read("/d/scratch/ASTR5160/week2/struc.fits")

def plotradec(RA,DEC):
    """
    This function plots Right Assencion and Declination from
    an "objs" file.
    """
    plt.plot(RA,DEC, "bx")
    
    
"""
#Overplot just he RA and DEC of the objects with ext > 0.22
for i in range(len(objs["RA"])):
    if float(objs["EXTINCTION"][:,0]) > 0.22:
        plotradec(objs["RA"],objs["DEC"])
"""

#Post submission correction, for future coding via Adam:
ii = objs["EXTINCTION"][:,0] > 0.22
plotradec(objs["RA"][ii],objs["DEC"][ii])


