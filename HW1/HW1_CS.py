#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
ASTR 5160

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
from astropy.coordinates import EarthLocation
from astropy.io import ascii
from tqdm import tqdm


#CS Code takes as input a month, and reads in Q list, 
#and prints out which can be observed at lowest airmass from KPNO at 11pm MST
#CSNon table global variables
raq = []
deq = []

#CS Global variables that will eventually go in table
rat = []
dect = []
times = []
airmass = []

#CS File location
fl = "/d/scratch/ASTR5160/week4/HW1quasarfile.txt"


def RfileR(fstring):
    """
    CS This function takes a file location as a string
    and extracts RA.
    """
    objs = open(fstring,"r")
    for o in objs:
        ra = str(o[0:2])+"h"+str(o[2:4])+"m"+str(o[4:9])+"s"
        raq.append(ra)
    print("RA extracted from", fl)
    return raq

def RfileD(fstring):
    """
    CS This function takes a file location as a string
    and extracts DEC.
    """
    objs = open(fstring,"r")
    for o in objs:
        dec = str(o[9:12])+"d"+str(o[12:14])+"m"+str(o[14:18])+"s"
        deq.append(dec)
    print("DEC extracted from", fl)
    return deq



def cyear():
    cyear = Time.now()
    cyear.format = "byear"
    cyd = str(cyear.decimalyear)
    return cyd[0:4]


#CS how many days are in the given month
d31 = [1,3,5,7,8,10,12]
d30 = [4,6,9,11]

def getNdays(mn):
    """
    CS This function gets the number of days in a month,
    for the current year, accounting for leap years
    """
    #CS get the current year
    cyd = cyear()
    print("The current year is:",cyd[0:4])

    days = 0
    #CS get nubmer of days
    for i in range(len(d31)):
        if mn == d31[i]:
            days = 31
    for i in range(len(d30)):
        if mn == d30[i]:
            days = 30
    if int(cyd[0:4])%4 != 0:
        if mn == 2:
            days = 2
    else:
        if mn == 2:
            days = 29
            print("This year is a leap year")
    print("Number of days in",mn,"is:",days,"days")
    return days


def getTimes(mn):
    """
    This function gets a list of times for each day in a given month,
    at 11pm MST and converts them to UTC
    """
    
    #CS Note that MST=UTC-7hrs
    daysnum = getNdays(mn)
    days = np.arange(2,(daysnum+1),1)
    #CS note that "6 am on the first" is really 11pm MST on the last day
    #of the previous month
    cyd = cyear()
    #CS compute 11PM MST on each day of the month
    for i in range(len(days)):
        """
        CS The MST time of 11pm on the last day of the month is 6am on the
        first day of the next month.
        """
        if days[i] == (daysnum):
            mnew = mn+1
            if mnew == 13:
                mnew = 1
                times.append(str(int(cyd[0:4])+1)+"-"+str(mnew)+"-"+str(days[i])+"T06:00:00")
            else:
                times.append(str(cyd[0:4])+"-"+str(mnew)+"-01T06:00:00")

        if 1 <= days[i] <= 9:
            times.append(str(cyd[0:4])+"-"+str(mn)+"-0"+str(days[i])+"T06:00:00")
        if 10 <= days[i] <= (daysnum+1):
            times.append(str(cyd[0:4])+"-"+str(mn)+"-"+str(days[i])+"T06:00:00")
    #print("11pm MST in UTC for each day:",times)
    return times


def calAmass(ra,dec,time):
    """
    CS this function calculates minimum airmass at KPNO
    of a given object at a specified time
    """
    #CS get the skycoords of the input object
    coords = SkyCoord(ra,dec,frame = 'icrs')
    KPNO = EarthLocation.of_site("kpno")
    
    #CS Convert the input skycoords to AltAz, where Coords is array like
    for i in tqdm(range(len(time))):
        aa = AltAz(location=KPNO,obstime=Time(time[i]))
        caa = coords.transform_to(aa)
        maltns = max(caa.alt)
        malt = str(maltns)
        print("max alt:",malt)
        maltd = float(malt[0:2])+float(malt[3:5])/60+float(malt[6:8])/60
        AM = (1/math.cos(float(maltd)*180/(6.28)))
        airmass.append(AM)
        for k in range(len(caa.alt)):
            if maltns == caa.alt[k]:
                print("caa corresponding to max:",caa)
    #CS now compute air mass
    
    
def WTab():
    """
    CS Writes an ascii table in the desired format, with coordinates,
    date, and airmass for a given month. 
    """
    data = Table()
    data["Date"] = times
    #data["Quasar Coordinates (hms.ss dms)"] = SkyCoord(raq,deq,frame='icrs')
    #data["RA (deg)"] = raq
    #data["Dec (deg)"] = deq
    data["Airmass"] = airmass
    ascii.write(data, 'Q_Out_Table.dat', overwrite=True)
    
if __name__ == "__main__": 
    #Rfile(fl)
    month = input("Input a desired month (1-12):")
    #getTimes(month)
    calAmass(RfileR(fl),RfileD(fl),getTimes(month))
    WTab()
    
    


