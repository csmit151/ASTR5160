# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 17:15:03 2025

@author: UW-User
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy import units as u

#Chase L. Smith
#Survey Observations

"""
Convert some RA (in the format deg,',") to deg and
DEC (in the format hms) to deg
"""
c = SkyCoord('02h37m23s','-37d17m32s', frame = 'icrs')
print("Converted Cooridantes are:",c,"Deg")

#We can check this against Adam's equations

#RA
RAhours = 2
RAmin = 37
RAsec = 23
#DEC
DECdeg = -37
Decmin = 17
Decsec = 32

def RA(h,m,s):
    print("Degrees RA:",15*(float(h)+float(m)/60+float(s)/3600))
def DEC(d,m,s):
    if float(d) < 0:
        print("Degrees DEC:",float(d)-float(m)/60-float(s)/3600)
    else:
        print("Degrees DEC:",float(d)+float(m)/60+float(s)/3600)


RA(RAhours,RAmin,RAsec)
DEC(DECdeg,Decmin,Decsec)

"""
Output:
    
Skycoord:
Converted Cooridantes are: <SkyCoord (ICRS): (ra, dec) in deg
    (39.34583333, -37.29222222)> Deg
    
Adam's eq:
Degrees RA: 39.34583333333334
Degrees DEC: -37.29222222222222

Tada!
"""

#print(Time.now())

#Todays Date in both JD and MJD
tjd = Time.now()
tjd.format = 'jd'
tmjd = Time.now()
tmjd.format = 'mjd'

print("Today in JD:", tjd)
print("Today in MJD", tmjd)

#And they are indeed diffrent by 2400000.5

#Use np.arange and Time.now to list some days near today
dfn = np.arange(0,10,1) #Days from now
for i in range(len(dfn)):
    t1 = Time.now()
    print(t1+dfn[i])