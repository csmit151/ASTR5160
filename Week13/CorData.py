#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fitting A Line: X^2 with correlated data

@author: csmit151
"""

# CS Import relevant packages
import numpy as np
import matplotlib.pyplot as plt

# CS Impot line.data
line_data = np.loadtxt("/d/scratch/ASTR5160/week13/line.data", unpack=True)

# CS Use np.cov to find the covariance matrix
co_var = np.cov(line_data)
print(co_var.shape)

# CS Resulting matrix is 10x10,
# CS which makes sense as there are
# CS 10 x bins, ie, sigma values run 1-10 by 1-10
# CS We want to know how data in one bin is correlated
# CS with data in another bin.
# CS You could take a million measurments in each bin and
# CS its not going to change how the underlying bins are correlated
# CS just the number of data points you have to measure the correlation with

var = []
# CS Calculate the variance and print diagonal covariance
for i in range(0,10,1):
    print("var:",np.var(line_data[i],ddof=1))
    print("co_var:",co_var[i][i])
    var.append(np.var(line_data[i]))


# CS Confirmed that the diagonal of the covar
# CS is the variance

# CS Determine the values and locations fo the most anti-correlated
# CS and the most correlated columns of data, ignoring the diagonal

max_cor = []
min_cor = []
co_var_nd = [] # CS Covariance no diagonal
col_row = [] # CS Column and row of the 

for i in range(0,10,1):
    for j in range(0,10,1):
        if i != j:
            co_var_nd.append(co_var[i][j])
            col_row.append((str(i)+" "+str(j)))

max_cor_ind = np.argmax(co_var_nd)
min_cor_ind = np.argmin(co_var_nd)


print("Maximum covar:",np.max(co_var))
print("Maximum covar ignoring diagonal:",np.max(co_var_nd))
print("Data columns that are most correlated are:", col_row[max_cor_ind])
print("Minimum covar (most anti-correlated):",np.min(co_var_nd))
print("Minimum covar ignoring diagonal:",np.min(co_var))









