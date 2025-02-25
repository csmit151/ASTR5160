#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
ASTR 5160

Chase L. Smith

HW 1
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
import re
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
timecor = [] #CS Only times that are positive airmass
airmass = []
radecout = []
raout = []
decout = []

#CS File location
fl = "/d/scratch/ASTR5160/week4/HW1quasarfile.txt"

if __name__ == "__main__":
    def RfileR(fstring):
        """
        CS This function takes a file location as a string
        and extracts RA.
        """
        objs = open(fstring,"r")
        #CS Go through each object in the text file
        for o in objs:
            #CS Extract the first two numbers as RA hours, the next two as min,
            #and the next 5 as seconds
            ra = str(o[0:2])+"h"+str(o[2:4])+"m"+str(o[4:9])+"s"
            #CS Put these into an array to use later
            raq.append(ra)
        #print("RA extracted from", fl)
        return raq

if __name__ == "__main__":
    def RfileD(fstring):
        """
        CS This function takes a file location as a string
        and extracts DEC.
        """
        objs = open(fstring,"r")
        for o in objs:
            #CS Extract the next three as +hh, the next two as min,
            #and the last 4 as seconds
            dec = str(o[9:12])+"d"+str(o[12:14])+"m"+str(o[14:18])+"s"
            #CS Put these into an array to use later
            deq.append(dec)
        #print("DEC extracted from", fl)
        return deq


if __name__ == "__main__": 
    def cyear():
        """
        CS This function returns the current year
        """
        #CS Get the current time
        cyear = Time.now()
        cyear.format = "byear"
        #CS Extract the current year of from the current time
        cyd = str(cyear.decimalyear)
        return cyd[0:4]


#CS how many days are in the given month
d31 = [1,3,5,7,8,10,12]
d30 = [4,6,9,11]

if __name__ == "__main__": 
    def getNdays(mni):
        """
        CS This function gets the number of days in a month,
        for the current year, accounting for leap years
        """
        #CS get the current year
        cyd = cyear()
        print("The current year is:",cyd[0:4])
    
        days = 0
        #CS get nubmer of days and
        #check if this month as 31 or 30 days,or is Feb
        mn = int(mni)
        
        if days == 0:
            for i in range(len(d31)):
                if mn == d31[i]:
                    days = 31
            for i in range(len(d30)):
                if mn == d30[i]:
                    days = 30
            #CS If the year is an even multiple of 4, then
            #the current year must be a leap year
            if int(cyd[0:4])%4 != 0:
                if mn == 2:
                    days = 28
            else:
                #Otherwise use the normal num of days for feb
                if mn == 2:
                    days = 29
                    print("This year is a leap year")
            print("Number of days in input month,",mn,", is:",days,"days")
        return days

if __name__ == "__main__": 
    def getTimes(mn1):
        """
        This function gets a list of times for each day in a given month,
        at 11pm MST and converts them to UTC
        """
        
        #CS Note that MST=UTC-7hrs
        daysnum = getNdays(mn1)
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
                #CS Go to the next month if it is the last MST day of the month
                mnew = int(mn1)+1
                if mnew == 13:
                    #CS Go to Jan of the next year if 11PM MST is on 12-31
                    mnew = 1
                    times.append(str(int(cyd[0:4])+1)+"-"+str(mnew)+"-"+str(days[i])+"T06:00:00")
                else:
                    times.append(str(cyd[0:4])+"-"+str(mnew)+"-01T06:00:00")
            #CS append days to times, adding appropriate num of zeros
            if 1 <= days[i] <= 9:
                times.append(str(cyd[0:4])+"-"+str(mn1)+"-0"+str(days[i])+"T06:00:00")
            if 10 <= days[i] <= (daysnum+1):
                times.append(str(cyd[0:4])+"-"+str(mn1)+"-"+str(days[i])+"T06:00:00")
        #print("11pm MST in UTC for each day:",times)
        return times

if __name__ == "__main__": 
    def calAmass(ra,dec,time):
        """
        CS this function calculates minimum airmass at KPNO
        of a given object at a specified time
        """
        #CS get the skycoords of the input object
        coords = SkyCoord(ra,dec,frame = 'icrs')
        KPNO = EarthLocation.of_site("kpno")
        
        #CS For each time convert the input skycoords to AltAz, where Coords is array like
        for i in range(len(time)):
            #CS Alt and Az of an object at KP at each time
            aa = AltAz(location=KPNO,obstime=Time(time[i]))
            caa = coords.transform_to(aa)
            #CS Find the maximum altiude object (minimum airmass)
            #CS Note that am will be negative for z > 90
            alts = str(caa.alt)
            maltns = max(caa.alt)
            malt = str(maltns)
            #CS Convert to degrees
            maltd = float(malt[0:2])+float(malt[3:5])/60+float(malt[6:8])/60
            #Calculate Airmass
            AM = (1/math.cos(float(maltd)*180/(3.14)))
            #CS we only want airmass greater than zero
            if AM > 0:
                airmass.append(AM)
                #CS Record the time and this objects coordinates
                timecor.append(time[i])
                for j in range(len(caa)):
                    if str(caa.alt[j]) == malt:
                        radecout.append(str(coords[j].ra)+", "+str(coords[j].dec))
                        raout.append(str(coords[j].ra.degree))
                        decout.append(str(coords[j].dec.degree))
            
            if AM <= 0:
                #CS Remove negative airmasses
                #CS Only running this when AM < 0 saves us from needlessly
                #calculating airmasses for every singel object, everysingle day
                AMall = [] #CS An AM list w/o negatives removed
                AMsort = [] #CS A list for sorting non negative AM
                #CS Calculate all airmasses
                
                for j in range(len(caa)):
                    #CS Convert to degree and calculate Airmass
                    alts2 = str(caa[j].alt.degree)
                    AMall.append(1/math.cos(float(alts2)*180/(3.14)))
                #CS Find positive AMs only
                for l in range(len(AMall)):
                    if AMall[l] > 0:
                        AMsort.append(AMall[l])
                #CS find minimum, while keeping track of index
                Amin = 10 #CS A very high airmass to start sorting with
                intAmin = 0 #CS the ith AM value that is minimized, so we can get ra and dec later
                for k in range(len(AMsort)):
                    if AMsort[k] <= Amin:
                        Amin == AMsort[k]
                        intAmin == k
                #CS Write our found min airmass and its corresponding coords
                airmass.append(Amin)
                timecor.append(time[i])
                for j in range(len(caa)):
                    if str(caa.alt[j]) == str(caa.alt[intAmin]):
                        radecout.append(str(coords[j].ra)+", "+str(coords[j].dec))
                        raout.append(str(coords[j].ra.degree))
                        decout.append(str(coords[j].dec.degree))

            
        
if __name__ == "__main__": 
    def WTab():
        """
        CS Writes an ascii table in the desired format, with coordinates,
        date, and airmass for a given month. 
        """
        data = Table()
        data["Date"] = timecor
        data["Quasar Coordinates (hms.ss dms)"] = radecout
        data["RA (deg)"] = raout
        data["Dec (deg)"] = decout
        data["Airmass"] = airmass
        ascii.write(data, 'Q_Out_Table.dat', overwrite=True)
    
if __name__ == "__main__":
    #CS Take an input month from the comand line
    month = input("Input a desired month (1-12):")
    #CS Run calculate airmass script
    print("Running...")
    calAmass(RfileR(fl),RfileD(fl),getTimes(month))
    print("Writing...")
    #CS Write to data file
    WTab()
    
    


