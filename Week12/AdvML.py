#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Adventures in Machine Learning

@author: csmit151
"""

import numpy as np
from sklearn.datasets import load_iris
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, join
import os
from tqdm import tqdm
from glob import glob
from Week8.CrossMatch import decode_sweep_name
from Week8.CrossMatch import is_in_box
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import fits
import random

iris = load_iris()
#print(iris.DESCR)



#print("First few rows of the data:")
colnames = "sepal_length sepal_width petal_length petal_width" 
#print(colnames)
#print(iris.data[:3])
#print("---------------------------------")
#print("Column statistics:")
#for i, j in enumerate(iris.data[0]):
#      print("column: {}, mean: {:.2f}".format(i, np.mean(iris.data[..., i])))
#print("---------------------------------")
#print("Total number of rows: {}".format(len(iris.data)))
#print("Target classes (which type of iris): {}".format(iris.target))
#print("Target class names {}".format({i:j for i, j in enumerate(iris.target_names)}))


import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 1, figsize=(8,6))
for i in range(3):
    target_class = iris.target == i
    ax.scatter(iris.data[target_class, 0], iris.data[target_class, 1], s=90, label=iris.target_names[i])
    ax.tick_params(labelsize=14)
    ax.set_xlabel("Sepal length (cm)", size=14)
    ax.set_ylabel("Sepal width (cm)", size=14)
    ax.legend(prop={'size': 14})
    
#plt.show()

from sklearn import neighbors
# ADM use the k-nearest neighbors algorithm (k-NN)
# ADM with distances only to the nearest neighbor.
knn = neighbors.KNeighborsClassifier(n_neighbors=1)
knn.fit(iris.data, iris.target)
#print(colnames)
mock_data = [5, 4, 1, 0]
#print(knn.predict([mock_data]), iris.target_names[knn.predict([mock_data])])
mock_data = [6, 3, 4, 1]
#print(knn.predict([mock_data]), iris.target_names[knn.predict([mock_data])])

# ADM let's map out the entire sepal_length/sepal_width space!
# ADM (I'm restricting to just sepal_length and sepal_width as
# ADM it's easier to picture a 2-D space than a 4-D space).
n = 100000
mock_data = []
# ADM normally I don't condone appending to empty lists, but here
# ADM I want to explicitly illustrate which columns I'm working
# ADM on. This won't be a slow append, as it's only two columns.
for i in range(2):
    #print("working on column: {}".format(colnames.split()[i]))
    col_min = np.min(iris.data[..., i])
    col_max = np.max(iris.data[..., i])
    # ADM generate random points in the space corresponding to the
    # ADM iris measurement of interest.
    mock_meas = np.random.random(n)*(col_max - col_min) + col_min
    mock_data.append(mock_meas)
# ADM we now have a list of n*2 measurements, 
# ADM but we want an array of 2 columns and n rows.
mock_data = np.reshape(mock_data, (2, n)).T
mock_data





# ADM classify using the k-NN "black box"
# ADM trained on the real-world iris data.
# ADM but only use 2 columns as a simple illustration.
knn = neighbors.KNeighborsClassifier(n_neighbors=1)
knn.fit(iris.data[..., :2], iris.target)
mock_target_class = knn.predict(mock_data)
# ADM again, this for loop isn't strictly necessary, but
# ADM it's a clear way to print the information to screen.
#for i in range(10):
#    print(mock_data[i], mock_target_class[i], iris.target_names[mock_target_class[i]])


# CS Approximately what precentage of thses 100,000 irises
# CS are classifed as virgincia by the k-NN algorithim?

ii = (mock_target_class == 1)
#print("% virginica:")
#print(float(len(mock_data[ii])/len(mock_data))*100,"%")

"""
CS Results
-------------------
% virginica:
21.869 %
"""

"""

# CS Copy from last class
# CS Find the closest object in the sweep files to 
ra = 188.53667
dec = 21.04572

objs = Table(names=("RA", "DEC"))
objs.add_row((ra,dec))

fns = glob("/d/scratch/ASTR5160/data/legacysurvey/dr9/south/sweep/9.0/*fits")

qs_to_match = Table.read("/d/scratch/ASTR5160/week10/qsos-ra180-dec30-rad3.fits")
"""

"""
# CS Loop through all of the sweep files and use is_in_box to
# CS  get the one that contains our object6
for i in range(len(fns)):
    ramin, ramax, decmin, decmax = decode_sweep_name(fns[i])
    radecbox = [ramin, ramax, decmin, decmax]
    ii = is_in_box(objs,radecbox)
    if np.any(ii):
        print(fns[i])
"""

"""
CS Our sweep file is:
    /d/scratch/ASTR5160/data/legacysurvey/dr9/south/sweep/9.0/sweep-180p020-190p025.fits
"""

#sweep1 = Table.read("/d/scratch/ASTR5160/data/legacysurvey/dr9/south/sweep/9.0/sweep-180p020-190p025.fits")

#ii = (sweep1["RA"] > (ra-0.0001)) & (sweep1["RA"] < (ra+0.0001)) & (sweep1["DEC"] > (dec-0.0001)) & (sweep1["DEC"] < (dec+0.0001))

#smallsweep = sweep1[ii]


# CS Retrive every point source in the sweeps that are within 3 deg of
# CS ra, dec = 180, 30

# CS Star by finding which sweep files we need
# CS Loop through all of the sweep files and use is_in_box to
# CS  get the one is in our box

"""
boxobj = Table(names=("RA", "DEC")) # Make a box that encloses 3 deg in all directions
boxobj.add_row((183,33))
boxobj.add_row((177,33))
boxobj.add_row((177,27))
boxobj.add_row((183,27))




sweeplist = []

for i in range(len(fns)):
    ramin, ramax, decmin, decmax = decode_sweep_name(fns[i])
    radecbox = [ramin, ramax, decmin, decmax]
    ii = is_in_box(boxobj,radecbox) # CS Find the sweep files we need to get all of points in the sweeps that are within 3deg of the location
    if np.any(ii):
        #print(fns[i])
        sweeplist.append(str(fns[i]))

# CS Get every point source with TYPE == PSF and r mag < 20
sw1 = []
sw2 = []
sw3 = []
sw4 = []

swtotlen = 0

for i in tqdm(range(len(sweeplist))):
    ra1, dec1 = [180.]*u.degree,[30.]*u.degree
    c1 = SkyCoord(ra1,dec1,frame='icrs')
    if i == 0:
        sweep = Table.read(sweeplist[i])
        ra2 = np.array(sweep["RA"])*u.degree
        dec2 = np.array(sweep["DEC"])*u.degree
        c2 = SkyCoord(ra2,dec2,frame='icrs')
        kk = (sweep["FLUX_R"] > 10) & (c1.separation(c2) < 3*u.degree)# CS Corresponds to r < 20
        swa = sweep[kk]
        ii = (swa["TYPE"] == "PSF")
        sw1 = swa[ii]
        print(len(sw1["RA"]))
        swtotlen = swtotlen + len(sw1["RA"])
    if i == 1:
        sweep = Table.read(sweeplist[i])
        ra2 = np.array(sweep["RA"])*u.degree
        dec2 = np.array(sweep["DEC"])*u.degree
        c2 = SkyCoord(ra2,dec2,frame='icrs')
        kk = (sweep["FLUX_R"] > 10) & (c1.separation(c2) < 3*u.degree)# CS Corresponds to r < 20
        swa = sweep[kk]
        ii = (swa["TYPE"] == "PSF")
        sw2 = swa[ii]
        print(len(sw2["RA"]))
        swtotlen = swtotlen + len(sw2["RA"])
    if i == 2:
        sweep = Table.read(sweeplist[i])
        ra2 = np.array(sweep["RA"])*u.degree
        dec2 = np.array(sweep["DEC"])*u.degree
        c2 = SkyCoord(ra2,dec2,frame='icrs')
        kk = (sweep["FLUX_R"] > 10) & (c1.separation(c2) < 3*u.degree)# CS Corresponds to r < 20
        swa = sweep[kk]
        ii = (swa["TYPE"] == "PSF")
        sw3 = swa[ii]
        print(len(sw3["RA"]))
        swtotlen = swtotlen + len(sw3["RA"])
    if i == 3:
        sweep = Table.read(sweeplist[i])
        ra2 = np.array(sweep["RA"])*u.degree
        dec2 = np.array(sweep["DEC"])*u.degree
        c2 = SkyCoord(ra2,dec2,frame='icrs')
        kk = (sweep["FLUX_R"] > 10) & (c1.separation(c2) < 3*u.degree)# CS Corresponds to r < 20
        swa = sweep[kk]
        ii = (swa["TYPE"] == "PSF")
        sw4 = swa[ii]
        print(len(sw4["RA"]))
        swtotlen = swtotlen + len(sw4["RA"])


print(swtotlen)

# CS Coordinate match psfobjs to the qsosfile

psfobjs = np.concatenate((sw1,sw2,sw3,sw4), axis=0)

qra1 = np.array(qs_to_match["RA"])*u.degree
qdec1 = np.array(qs_to_match["DEC"])*u.degree
c1 = SkyCoord(qra1,qdec1)

pra2 = np.array(psfobjs["RA"])*u.degree
pdec2 = np.array(psfobjs["DEC"])*u.degree
c2 = SkyCoord(pra2,pdec2,frame='icrs')

id1, id2, d2, d3 = c2.search_around_sky(c1, (1/3600)*u.deg) # CS 1" mr

"""


"""
# CS id1 is 275 r<20 quasars

# CS write a function that uses the k-NN algorithm to
# CS classify quasars based on g-z and r-W1 colors

print("Now writing table")


with open("id1.txt", "w") as file:
    for i in range(len(id1)):
        file.write(str(id1[i]))
        file.write("\n")

with open("id2.txt", "w") as file:
    for i in range(len(id2)):
        file.write(str(id2[i]))
        file.write("\n")

print("Table written")

"""



"""
# CS Objects we know are qusars

match_psf = psfobjs[id1]  
print(len(match_psf["RA"]))


# CS write out relevant data for ease of input, no problem as only 275 objects
with open("fluxR.txt", "w") as file:
    for i in range(len(match_psf["FLUX_R"])):
        file.write(str(match_psf["FLUX_R"][i]))
        file.write("\n")

with open("fluxW1.txt", "w") as file:
    for i in range(len(match_psf["FLUX_W1"])):
        file.write(str(match_psf["FLUX_W1"][i]))
        file.write("\n")

with open("fluxZ.txt", "w") as file:
    for i in range(len(match_psf["FLUX_Z"])):
        file.write(str(match_psf["FLUX_Z"][i]))
        file.write("\n")
        
with open("fluxG.txt", "w") as file:
    for i in range(len(match_psf["FLUX_G"])):
        file.write(str(match_psf["FLUX_G"][i]))
        file.write("\n")

"""

#print(len(pdata["RA"]))

fluxr = open('fluxR.txt') # CS read in id1
fr1 = fluxr.read()
fr = fr1.splitlines()
fluxr.close()

fluxw1 = open('fluxW1.txt') # CS read in id1
fw11 = fluxw1.read()
fw1 = fw11.splitlines()
fluxw1.close()

fluxg = open('fluxG.txt') # CS read in id1
fg1 = fluxg.read()
fg = fg1.splitlines()
fluxg.close()

fluxz = open('fluxZ.txt') # CS read in id1
fz1 = fluxz.read()
fz = fz1.splitlines()
fluxz.close()

g_z = []
r_w1 = []

for i in range(len(fr)):
    g_z.append(float(fg[i])-float(fz[i]))
    r_w1.append(float(fr[i])-float(fw1[i]))



"""
# CS Use a random subsample of psfobjs as "stars"
rss = []

for i in range(0,274,1):
    randint = random.randint(0, (len(psfobjs)-1))
    rss.append(randint)

rand_psf = psfobjs[rss]



with open("rand_fluxR.txt", "w") as file:
    for i in range(len(rand_psf["FLUX_R"])):
        file.write(str(rand_psf["FLUX_R"][i]))
        file.write("\n")

with open("rand_fluxW1.txt", "w") as file:
    for i in range(len(rand_psf["FLUX_W1"])):
        file.write(str(rand_psf["FLUX_W1"][i]))
        file.write("\n")

with open("randfluxZ.txt", "w") as file:
    for i in range(len(rand_psf["FLUX_Z"])):
        file.write(str(rand_psf["FLUX_Z"][i]))
        file.write("\n")
        
with open("rand_fluxG.txt", "w") as file:
    for i in range(len(rand_psf["FLUX_G"])):
        file.write(str(rand_psf["FLUX_G"][i]))
        file.write("\n")

"""

#print(len(pdata["RA"]))

rand_fluxr = open('rand_fluxR.txt') # CS read in id1
rand_fr1 = rand_fluxr.read()
rand_fr = rand_fr1.splitlines()
rand_fluxr.close()

rand_fluxw1 = open('rand_fluxW1.txt') # CS read in id1
rand_fw11 = rand_fluxw1.read()
rand_fw1 = rand_fw11.splitlines()
rand_fluxw1.close()

rand_fluxg = open('rand_fluxG.txt') # CS read in id1
rand_fg1 = rand_fluxg.read()
rand_fg = rand_fg1.splitlines()
rand_fluxg.close()

rand_fluxz = open('randfluxZ.txt') # CS read in id1
rand_fz1 = rand_fluxz.read()
rand_fz = rand_fz1.splitlines()
rand_fluxz.close()

rand_g_z = []
rand_r_w1 = []

for i in range(len(rand_fr)):
    rand_g_z.append(float(rand_fg[i])-float(rand_fz[i]))
    rand_r_w1.append(float(rand_fr[i])-float(rand_fw1[i]))


"""
plt.scatter(g_z,r_w1, color='b')
plt.scatter(rand_g_z,rand_r_w1, color='r')

plt.xlabel("g_z")
plt.ylabel("r_w1")
plt.show()
"""


    
"""
# ADM let's plot the sepal_length, sepal_width space of the k-NN classifier.
fig, ax = plt.subplots(1, 1, figsize=(8,6))
for i in range(3):
    target_class = mock_target_class == i
    ax.scatter(mock_data[target_class, 0], mock_data[target_class, 1], s=10, label=iris.target_names[i])
    ax.tick_params(labelsize=14)
    ax.set_xlabel("Sepal length (cm)", size=14)
    ax.set_ylabel("Sepal width (cm)", size=14)
    ax.legend(prop={'size': 14})
    plt.show()
"""  
    
# CS Now to astronomy:
    

qsos = np.array([
        [18.744722,18.394268,18.309591,18.184225,18.261715],[19.688782,19.054688,19.176489,19.02421,18.696566],[19.563297,19.471746,19.38908,19.067038,19.207346],[19.029963,18.972807,18.958265,18.823439,18.692583],[18.149904,18.051216,18.112373,18.269794,17.576611],[19.233286,19.249971,19.067638,19.106289,18.923378],[19.342098,19.114586,19.048199,18.923874,19.003283],[19.572979,19.539234,19.593983,19.58787,19.839739],[19.451731,19.517733,19.291042,18.968428,19.06872],[19.305523,18.724361,18.139717,17.776787,17.693134],[19.358782,19.263777,19.133242,18.959019,18.794546],[19.572119,19.517012,19.309811,19.112907,19.153662],[19.551817,19.186308,18.852896,18.728289,18.694551],[18.16538,18.106567,18.094177,17.904747,17.910076],[19.453842,19.284033,19.102739,18.83543,18.96608],[19.854441,19.243029,18.528109,18.219172,18.114115],[20.241364,20.597746,20.206196,20.01281,19.859661],[19.862375,19.731724,19.702806,19.616152,19.490868],[18.845604,18.700039,18.667009,18.467548,18.478443],[19.482214,19.203161,19.032568,18.737476,18.744768],[19.3901,19.198641,19.013081,19.064732,18.865631],[20.719276,20.281059,19.662498,19.624067,19.32023],[19.56996,18.999855,18.937254,18.869497,18.590803],[17.922159,17.750494,17.517071,17.523623,17.361856],[21.113659,19.708879,18.911526,18.615511,18.395302],[22.712126,20.320133,19.86651,19.737288,19.70672],[19.800406,18.592144,17.817379,17.399628,17.060135],[22.647556,21.77388,20.735937,20.335043,19.775301],[20.881166,20.447176,20.348919,20.088438,20.099693],[19.402994,19.042908,18.94466,18.93589,18.843378],[18.903831,18.572012,18.352116,18.158442,18.000324],[19.914022,19.80702,19.525486,19.125696,18.892767],[19.332512,19.08815,18.854038,18.616337,18.555649],[19.259775,18.842321,18.864609,18.796721,18.658291],[18.393356,18.218657,18.018612,17.694494,17.581081],[18.694517,18.643589,18.597673,18.65523,18.452408],[20.336267,20.060621,19.405148,19.127205,18.881546],[18.658436,18.392712,18.137266,18.08651,18.109537],[19.511658,19.339643,19.001104,19.020147,19.032293],[19.610483,19.518547,19.175377,19.084261,18.924337],[20.707575,20.182711,19.369616,19.060297,18.839787],[19.601027,18.958122,18.881411,18.8493,18.668055],[19.532183,19.280308,19.297268,19.07955,19.029108],[18.781235,18.52058,18.253962,18.253763,18.105644],[20.283949,20.171352,20.088402,19.87583,19.820202],[19.820608,19.679596,19.588165,19.389709,19.168039],[22.094894,21.437128,20.592884,20.062315,19.739923],[19.526564,19.272533,19.07749,19.07902,18.948053],
        [18.902395,18.322453,17.777796,17.338715,17.216942],[20.018803,19.598724,19.489754,19.370985,19.299452],[21.604231,19.705097,19.355795,18.941904,18.658649],[21.492113,20.53562,20.504961,20.145609,19.690453],[20.506151,19.143932,18.914755,18.733116,18.68424],[19.779257,19.587122,19.346659,19.172829,19.228399],[19.080269,18.923189,18.714451,18.805676,18.715391],[23.893023,20.949961,20.108637,19.602491,19.402994],[19.475689,19.497339,19.21969,19.232433,19.362291],[20.722937,20.438725,20.412531,20.209305,20.110937],[20.378387,19.978376,19.957357,19.798298,19.801628],[19.394093,19.159006,19.039072,18.823759,18.833567],[18.6556,18.559763,18.525602,18.279596,18.298143],[19.82836,19.361555,19.132957,19.00437,18.755081],[18.314793,17.968035,17.770653,17.55607,17.568619],[20.738749,20.25322,20.429714,20.134119,19.80547],[18.92709,18.789314,18.328901,17.930676,17.71903],[20.765722,20.425152,20.261223,20.032152,19.987698],[21.898718,20.366604,19.190598,18.685432,18.564209],[19.209698,18.843922,18.082823,17.598646,17.454884],[20.571241,19.715168,19.346649,19.035482,18.673683],[19.700674,19.372719,19.263952,18.997389,18.928925],[19.755407,19.536802,19.430946,19.222158,19.370178],[19.109722,18.664597,18.477087,18.20159,18.220667],[22.367666,20.430206,19.458319,19.31806,19.301077],[20.666672,20.259563,19.443361,19.033512,18.698406],[18.640318,18.900045,18.630014,18.609161,18.669069],[19.134163,19.020254,18.747284,18.952831,19.055578],[21.476528,19.508688,19.135096,19.140932,19.117966],[19.61648,19.265966,19.163887,19.116779,19.071653],[18.611923,18.521744,18.484806,18.22971,18.105438],[19.663836,18.972357,18.550348,18.392038,18.211784],[19.677265,19.40823,19.089659,18.671545,18.642952],[19.203194,19.177202,19.019161,18.976236,19.129333],[19.157156,19.094862,19.058044,18.95915,18.794933],[17.064066,16.845049,16.747282,16.760361,16.699829],[19.263111,19.05188,18.779545,18.766621,18.74962],[21.415379,18.908577,18.490671,18.3186,18.210991],[19.178915,18.984844,18.734407,18.681,18.757259],[19.153744,18.989044,18.858898,18.750517,18.556971],[21.874535,19.867119,19.67396,19.572515,19.232058],[19.40521,19.081427,18.945744,18.801537,18.811073],[18.206612,17.959349,17.774895,17.571123,17.494965],[19.324575,19.241135,18.977043,18.90477,19.002703],[19.256691,18.664949,18.208923,17.746832,17.563005],[18.234364,17.920656,17.877764,17.82852,17.7251],[19.474617,19.370167,19.275578,19.164686,18.90004],[19.953865,19.784174,19.708263,19.449259,19.312752],[19.755163,19.178165,18.529041,18.113386,17.857853],[19.727612,19.195156,19.101727,19.117582,19.022566],[19.894958,19.867922,19.807981,19.763544,19.343872],[19.901382,19.600948,19.302818,19.175266,18.703302]
    ])



stars = np.array([
        [22.44946,19.97534,18.54526,17.33972,16.67254],[20.77815,19.4398,18.87869,18.71637,18.68406],[22.61502,20.0222,18.53809,17.37034,16.69467], [22.67518,20.3125,18.82124,17.89092,17.40948], [23.37306,20.26618,18.75942,17.93758,17.4218], [20.87663,19.24346,18.82649,18.61213,19.03254], [21.2547,19.44574,18.59049,18.26796,18.11473], [21.14091,19.84629,19.35966,19.11529,19.03632], [23.34694,20.00362,18.47734,17.81385,17.38007], [20.54105,19.3501,18.90781,18.71399,18.53207], [20.81281,19.75167,19.37374,19.15161,19.11015], [21.19628,19.91164,19.38225,19.06724,19.00649], [23.42917,20.47683,19.41293,18.92953,18.66926], [22.23457,19.4957,17.97546,17.33125,16.91896], [21.84833,20.16281,18.61554,17.70773,17.20546], [21.87615,19.51687,18.07097,17.28559,16.8108], [20.15546,19.30957,18.97326,18.88849,18.87114], [20.20057,19.31459,18.98439,18.86432,18.8286], [20.25109,19.29368,18.96605,18.85745,18.78543], [22.06602,19.48195,18.0576,17.27881,16.86323], [22.10527,19.50775,18.0136,17.24724,16.80717], [20.18234,19.29664,18.97465,18.86959,18.75895], [21.84393,19.52566,18.0703,17.25718,16.82474], [21.78772,19.53626,18.22131,16.77889,16.10115], [20.44573,19.01728,18.42748,18.23824,18.12197], [20.53322,19.02352,18.44298,18.22176,18.1349], [20.59585,19.01741,18.41915,18.24291,18.08544], [20.75146,19.64077,19.28507,19.11219,19.0484], [20.57132,19.30027,18.72266,18.55597,18.40891], [21.2337,19.77214,19.02792,18.81852,18.69435], [20.70082,19.25805,18.75136,18.56669,18.42665], [20.70589,19.8298,19.47613,19.40005,19.39293], [21.6109,19.78855,19.02922,18.81845,18.6857], [22.80736,20.32154,18.76478,16.92489,15.8764], [21.56886,19.43514,18.27148,17.80843,17.53949], [20.34764,19.34764,19.07848,18.94273,18.89366], [21.50406,19.62864,18.52377,17.95594,17.628], [21.37426,19.94507,19.23395,18.91821,18.75806], [20.9035,19.08836,18.30525,18.0126,17.84013], [21.71455,20.42793,19.44023,18.93189,18.62866], [20.86752,19.06702,18.29345,18.00766,17.82436], [21.98555,19.52187,18.51127,18.09584,17.88142], [21.67668,19.5505,18.50358,18.10969,17.87844],[22.53292,19.7276,18.48423,17.98095,17.70833],[20.8909,19.0599,18.28369,18.00476,17.87395],[22.35347,20.25845,18.95716,18.49136,17.97949],[21.76183,19.88117,18.75925,18.26558,17.91661],[22.23183,19.74561,18.488,17.97591,17.71163],[21.29818,19.79457,19.14499,18.88091,18.72307], [22.0105,19.70986,18.52199,17.97292,17.61687],
        [22.86553,19.94115,18.7831,18.23227,17.90711],[21.18186,19.79757,19.13424,18.85658,18.71214],[20.0074,19.05439,18.69278,18.60722,18.52868],[21.58379,20.33799,19.66424,19.2912,19.09761],[21.98114,20.2121,19.55524,19.32269,19.07644],[22.81192,19.89781,18.5065,17.95065,17.61704],[21.83415,20.23835,19.52526,19.28796,19.1398],[22.63741,20.26673,19.21909,18.85295,18.62377],[19.89656,19.03701,18.70751,18.5632,18.56942],[21.30959,19.37278,18.31021,17.86578,17.56223],[21.22793,19.36349,18.30933,17.81932,17.51374],[21.14737,19.33545,18.30335,17.81113,17.4971],[22.21014,20.27786,19.2433,18.83222,18.65571],[20.18058,19.31195,19.06964,18.94653,18.89464],[22.35099,20.28571,19.28649,18.80538,18.62735],[20.24579,19.3125,19.05523,18.95683,18.98032],[21.80952,20.44154,20.00894,19.91379,19.79607],[22.21284,19.61759,18.24859,17.67601,17.37054],[21.68078,20.48232,20.02647,19.88837,19.78827],[23.75356,20.39333,19.35249,18.94558,18.62885],[23.42696,20.37005,19.37519,18.97837,18.63843],[22.65762,19.64814,18.25478,17.70169,17.37972],[21.21484,20.05256,19.62379,19.53748,19.36399],[21.19142,20.05359,19.71335,19.56334,19.70693],[20.98901,20.05082,19.68979,19.58001,19.68044],[20.77375,20.00404,19.64923,19.53775,19.54221],[21.86685,19.9237,19.11284,18.88337,18.71919],[21.93681,19.90158,19.12132,18.86649,18.70179],[21.70023,20.09205,19.14815,18.7781,18.55863],[21.79366,20.03435,19.14712,18.70042,18.51328],[21.62113,19.99053,19.14881,18.94777,18.46128],[20.99038,19.95668,19.61966,19.47392,19.42933],[22.07217,19.46344,18.08623,17.07386,16.53181],[20.9025,19.97484,19.64865,19.49392,19.34935],[20.89883,19.97347,19.66507,19.40612,19.44094],[21.04365,19.97419,19.63968,19.50207,19.40869],[21.00765,19.95171,19.63981,19.48675,19.35925],[22.15741,19.48624,18.06823,17.09778,16.57183],[22.16418,19.77536,18.5058,18.01847,17.66681],[22.94029,19.56238,18.00614,16.93008,16.35485],[21.88229,20.49575,19.87854,19.59909,19.45779],[21.98588,20.47682,19.89064,19.59884,19.41466],[22.5875,19.57998,18.01291,16.93524,16.34619],[24.41736,19.71074,18.24074,17.58002,17.19191],[22.06751,19.64084,18.23813,17.5793,17.18487],[21.71258,19.76319,18.92926,18.4789,18.18599],[21.07738,19.50764,18.88156,18.57879,18.49265],[21.55381,19.24872,18.14714,17.67999,17.32854],[21.64327,19.25232,18.14592,17.67702,17.38654],[21.04674,19.48493,18.854,18.60168,18.56519]
    ])



data = np.concatenate([stars, qsos]) # CS Remember as a good example of how to use concatonate
data_class = np.concatenate([np.zeros(len(stars), dtype='i'), 
                            np.ones(len(qsos), dtype='i')])
target_names = np.array(["star", "QSO"])

# CS This makes a table of 0,1 corresponding to stars, and qsos respectivly


"""
# ADM let's plot some of the qso and star colors.
fig, ax = plt.subplots(1, 1, figsize=(8,6))
for i in range(2):
    target_class = data_class == i
    ax.scatter(data[target_class, 2], data[target_class, 3], s=30, label=target_names[i]) # CS Note quick label making thingy, no mpatches required
    ax.tick_params(labelsize=14)
    ax.set_xlabel("r-band magnitude", size=14)
    ax.set_ylabel("i-band magnitude", size=14)
    ax.legend(prop={'size': 14})
    plt.show()

"""
    
# ADM classify using the k-NN "black box"
# ADM trained on the real-world astro data.
knn = neighbors.KNeighborsClassifier(n_neighbors=1)
knn.fit(data, data_class)
qsolike = [19.03,18.97,18.95,18.82,18.69]
starlike = [20.00,19.05,18.69,18.61,18.53]
predicted_class = knn.predict([qsolike])
print(predicted_class, target_names[predicted_class])
predicted_class = knn.predict([starlike])
print(predicted_class, target_names[predicted_class])

# CS fix mistake that I only grabbed 274, not 275
rand_g_z.append(rand_g_z[0])
rand_r_w1.append(rand_r_w1[0])


# CS Data:
data_stars = []
data_qso = []
for i in range(len(rand_g_z)):
    data_stars.append([rand_g_z[i],rand_r_w1[i]])
    data_qso.append([g_z[i],r_w1[i]])


allfluxes = np.concatenate([np.array(data_stars),np.array(data_qso)]) # CS put all of the fluxes together
# CS Following ADM example:
aflux_class = np.concatenate([np.zeros(len(data_stars), dtype='i'), np.ones(len(data_qso), dtype='i')])
tnames = np.array(["star","QSO"]) # CS make sure matches our order



def knn_qso(aflux):
    knn = neighbors.KNeighborsClassifier(n_neighbors=1)
    knn.fit(aflux, aflux_class)
    for i in range(len(g_z)):
        thisarray = [g_z[i],rand_g_z[i]]
        print(thisarray)
        predicted_class = knn.predict([thisarray]) # CS QSO like
        print(predicted_class, tnames[predicted_class])



if __name__ == "__main__":
    knn_qso(allfluxes)
    print("% virginica:")
    print(float(len(mock_data[ii])/len(mock_data))*100,"%")
    
    
    
"""
Results
---------------
CS A few ids:
[1] ['QSO']
[-47.4999593, -88.10191400000001]
[1] ['QSO']
[-24.750729000000007, -87.53321]
[1] ['QSO']
[-96.316645, -9.931935000000001]
[1] ['QSO']
[-480.87037, -228.85047000000003]
[1] ['QSO']
[-60.222306, -174.80122899999998]
[1] ['QSO']
[-126.75318699999998, -163.67246999999998]
[0] ['star']
[-126.35480000000001, -102.964123]
[1] ['QSO']
[-4.525053, -63.04858]
[1] ['QSO']
[-73.653614, -674.8139000000001]
[0] ['star']

CS And the precentage virginica:
% virginica:
21.511 %

"""