#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fitting a Line

@author: csmit151
"""

# CS Import relevant packages
import numpy as np
import matplotlib.pyplot as plt


# CS Impot line.data
line_data = np.loadtxt("/d/scratch/ASTR5160/week13/line.data", unpack=True)

# CS Global variables
mean = []
var = []

# CS Calculate mean and variance
for i in range(0,10,1):
    mean.append(np.mean(line_data[i]))
    var.append(np.var(line_data[i]))

xbins = np.arange(0.5,10.5,1)

colors = ['r','b','g','c','m']*2 # CS Make uniform colors for easy viewing

# CS Scatter plot each point
fig, ax = plt.subplots(1, 1, figsize=(8,6))
for i in range(len(xbins)):
    ldat = line_data[i]
    for j in range(len(ldat)):
        ax.scatter(xbins[i],ldat[j],color=colors[i])
        ax.tick_params(labelsize=14)
        ax.set_xlabel("X", size=14)
        ax.set_ylabel("Y", size=14)

# CS Guess a range of m and b value that could fit the data
x_lin = np.linspace(0,10)
y = 3*x_lin+4.5
ax.plot(x_lin,y, color="red")

# CS Dosen't look too bad, so m = 2.5-3.5, b = 4-5
# CS Determine Chi**2 for a grid of m and b that is over this range

# CS Make a grid
m_grid = np.arange(2.5,3.5,0.01)
b_grid = np.arange(4,5,0.01)


# CS Global grid variables
Chi_grid = []
Chi_grid_m = []
Chi_grid_b = []


def Chi_grid_mb(mg,bg,xb,ldat,vr):
    for i in range(len(mg)):
        for j in range(len(bg)):
            # CS Calculate for each grid point
            Chi_i = []
            m = mg[i]
            for k in range(len(xb)):
                E = (float(m)*float(xb[k]))+bg[j] # CS Calculate expected for this bin and grid point
                O = ldat[k]
                for l in range(len(O)):
                    Chi_i.append((O[l]-E)**2/(vr[k])**2) # Calculate (O_i-E_i)^2/sig_i^2
            Chi_grid.append(np.sum(Chi_i)) # CS Summation
            Chi_grid_m.append(m) # CS Record the grid point this Chi^2 goes to
            Chi_grid_b.append(b_grid[j])
        

if __name__ == "__main__":
    print("Mean of observed y values for x1-x9:")
    print(mean)
    print("Variance of the y values for x1-x9")
    print(var)
    Chi_grid_mb(m_grid,b_grid,xbins,line_data,var)
    bfit = np.argmin(Chi_grid)
    m_best = Chi_grid_m[bfit]
    b_best = Chi_grid_b[bfit]
    print("Best fit m:",m_best)
    print("best fit b:",b_best)
    x_lin_b = np.linspace(0,10)
    y_b = m_best*x_lin+b_best
    ax.plot(x_lin_b,y_b, color = "green")
    plt.show()
