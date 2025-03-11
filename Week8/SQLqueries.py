#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Class 15 Tasks

Chase L. Smith
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table

if __name__ == "__main__":
    objs = Table.read("SDSS1.csv")
    i = (objs["g"]) <= 17
    k = (objs["g"] > 17) & (objs["g"] < 19)
    j = (objs["g"]) >= 19
    plt.scatter(objs["ra"][i],objs["dec"][i], marker = "o", s = 50)
    plt.scatter(objs["ra"][k],objs["dec"][k], marker = "o", s = 20, color = "red")
    plt.scatter(objs["ra"][j],objs["dec"][j], marker = "o", s = 10, color = "green")
    plt.show()


