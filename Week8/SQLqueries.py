#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Class 15 Tasks

Chase L. Smith
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table

def plotbymag(file):
    """Take an input file of ra, dec, and g mag and plot
    larger magnitudes as smaller circles in a scatter plot
    
    Parameters
    ----------
    file: :class: 'string'
        The location of the input file
    
    Returns
    -------
    A scatter plot    
    """
    
    objs = Table.read(str(file)) # CS read in the file
    i = (objs["g"]) <= 17 # CS seperate objects into three magnitude groups
    k = (objs["g"] > 17) & (objs["g"] < 19)
    j = (objs["g"]) >= 19
    # CS plot those groups
    plt.scatter(objs["ra"][i],objs["dec"][i], marker = "o", s = 50)
    plt.scatter(objs["ra"][k],objs["dec"][k], marker = "o", s = 20, color = "red")
    plt.scatter(objs["ra"][j],objs["dec"][j], marker = "o", s = 10, color = "green")
    plt.show()
    
if __name__ == "__main__":
    plotbymag("SDSS1.csv")


