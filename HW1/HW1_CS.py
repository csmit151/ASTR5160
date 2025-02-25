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
timecor = [] #Only times that are positive airmass
airmass = []

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
        print("RA extracted from", fl)
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
        print("DEC extracted from", fl)
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
        print("Days before loop:",days)
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
                    days = 2
                    print("This many days:",days)
            else:
                if mn == 2:
                    days = 29
                    print("This year is a leap year")
            print("Number of days in",mn,"is:",days,"days")
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
        
        #CS Convert the input skycoords to AltAz, where Coords is array like
        for i in tqdm(range(len(time))):
            #CS Alt and Az of an object at KP at each time
            aa = AltAz(location=KPNO,obstime=Time(time[i]))
            caa = coords.transform_to(aa)
            #CS Find the maximum altiude object (minimum airmass)
            #CS Note that am will be negative for z > 90
            alts = str(caa.alt)
            #CS Convert the altiude to a string, and check to see
            #if that z value is less than 90
            for j in range(2,len(alts),20):
                print("j:",j)
                print("J+1:",j+2)
                print("alts[(j):(j+2)]:",alts[(j):(j+2)])
                print("alts[(j-3):(j+5)]:",alts[(j-3):(j+5)])

                if alts[(j):(j+2)] != "..":
                    if alts[(j+1):(j+2)] == "s" or alts[(j+1):(j+2)] == "d":
                    #CS Skip bad coords?
                        if float(alts[(j+4):(j+5)]) < 90:
                            maltns = max(caa.alt)
                            malt = str(maltns)
                            #print("max alt:",malt)
                            #CS Convert to degrees
                            maltd = float(malt[0:2])+float(malt[3:5])/60+float(malt[6:8])/60
                            #Calculate Airmass
                            AM = (1/math.cos(float(maltd)*180/(3.14)))
                            airmass.append(AM)
                            #CS Remove this time so it dosen't get counted needlessly
                            timecor.append(time[i])
                    elif float(alts[(j):(j+2)]) != "s":
                        if float(alts[(j):(j+2)]) != "s]":
                            if float(alts[(j):(j+2)]) < 90:
                                maltns = max(caa.alt)
                                malt = str(maltns)
                                #print("max alt:",malt)
                                #CS Convert to degrees
                                maltd = float(malt[0:2])+float(malt[3:5])/60+float(malt[6:8])/60
                                #Calculate Airmass
                                AM = (1/math.cos(float(maltd)*180/(3.14)))
                                airmass.append(AM)
                                #CS Remove this time so it dosen't get counted needlessly
                                timecor.append(time[i])
                            
                            #for k in range(len(caa.alt)):
                            #    if maltns == caa.alt[k]:
                            #        print("caa corresponding to max:",caa)
        
if __name__ == "__main__": 
    def WTab():
        """
        CS Writes an ascii table in the desired format, with coordinates,
        date, and airmass for a given month. 
        """
        data = Table()
        data["Date"] = timecor
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
    
    


