#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ASTR 5160 Final Project

@author: csmit151
"""

# CS Import needed packages
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
import corner
import emcee
import matplotlib.patches as mpatches
from scipy.optimize import minimize


def emcee_line_fit(data_location):
    """ Uses a Bayseian framework to fit a linear
    model of the form y=mx+b to the data 
    
    Parameters
    ----------
    data_location: :class: 'string file-path'
        File-path to data file with four columns, x,s,sig,yerr
    Returns
    -------
    A corner plot of parameters' posterior probability distributions
    A plot of y(x) including error bars
    Prints to screen best fitting parameters and uncertanties
    based on the 16th, 50th, and 84th precentiles
    """
    # CS Import dataxy
    data = Table.read(str(data_location))
    x = data["x"]
    y = data["y"]
    yerr = data["yerr"] # CS estimated error (sigma) on the y data
    # CS Define a log likelihood function
    def log_likelihood(theta, x, y, yerr):
        """ Computes the log likelihood of input parameters contained in theta
        for data parameterized by x, y and yerr
        
        Parameters
        ----------
        theta: :class: 'array'
            Input function parameters (i.e. m and b for a line)
        x: :class: 'array'
            Data x values
        y: :class: 'array'
            Data y values
        yerr: :class: 'array'
            Estimated error on the y data
        Returns
        -------
        Log liklihood for given parameters and data
        """
        m, b = theta
        model = m * x + b # CS linear model
        sigma2 = yerr**2 + model**2
        return -0.5 * np.sum((y - model) ** 2 / sigma2 + np.log(sigma2))
    
    # CS Choose some inital values that make sense for our line
    m_true = -2
    b_true = 2
        
    np.random.seed(42)
    nll = lambda *args: -log_likelihood(*args)
    initial = np.array([m_true, b_true]) + 0.1 * np.random.randn(2)
    soln = minimize(nll, initial, args=(x, y, yerr))
    m_ml, b_ml = soln.x
    
    
    # CS Make our prior bounds realistic for our line
    def log_prior(theta):
        """ Computes a uniform prior for an input array of parameters
        
        Parameters
        ----------
        theta: :class: 'array'
            Input function parameters (i.e. m and b for a line)
            
        Returns
        -------
        0 if parameters are within prior range, and -inf otherwise
        """
        m, b = theta
        if -6 < m < 0 and -3 < b < 6.0:
            # CS A uniform prior is zero when parameters are withinoru specified range
            return 0.0
        return -np.inf
    
    def log_probability(theta, x, y, yerr):
        """ Combines log_prior and log_likelihood functions
        to comput log probability
        
        Parameters
        ----------
        theta: :class: 'array'
            Input function parameters (i.e. m and b for a line)
        x: :class: 'array'
            Data x values
        y: :class: 'array'
            Data y values
        yerr: :class: 'array'
            Estimated error on the y data
        Returns
        -------
        log_probability as log_prior + log_likelihood
        """
        lp = log_prior(theta)
        if not np.isfinite(lp):
            return -np.inf
        return lp + log_likelihood(theta, x, y, yerr)
    
    # CS Start the walkers near the maximum likelihood
    pos = soln.x + 1e-4 * np.random.randn(32, 2)
    nwalkers, ndim = pos.shape
    
    sampler = emcee.EnsembleSampler(
        nwalkers, ndim, log_probability, args=(x, y, yerr)
    )
    # CS 5000 walkers should be fine according to emcee tutorials, check burn in though
    sampler.run_mcmc(pos, 5000, progress=True);
    
    # CS plot the time series of the parameters in the chain
    fig, axes = plt.subplots(2, figsize=(10, 7), sharex=True)
    samples = sampler.get_chain()
    labels = ["m", "b"]
    for i in range(ndim):
        ax = axes[i]
        ax.plot(samples[:, :, i], "k", alpha=0.3)
        ax.set_xlim(0, len(samples))
        ax.set_ylabel(labels[i])
        ax.yaxis.set_label_coords(-0.1, 0.5)
    
    axes[-1].set_xlabel("step number");
    
    # CS Print steps required to burn in
    tau = sampler.get_autocorr_time()
    print("Burn in:",tau)
    
    # CS Burn in takes like 40 something steps so 100 should be fine to discard
    # CS Remove the burn in steps and thin by 15 steps
    flat_samples = sampler.get_chain(discard=100, thin=15, flat=True)
    print("flat list of samples shape:",flat_samples.shape)
    
    # CS Note that removing the burn in steps lets the emcee walkers actually have time
    # CS to move to a "good" location in our parameter space before we start "trusting" the values
    # CS that they give us. 
    
    fig = corner.corner(flat_samples, labels=labels)
        
    # CS Get fit params for easy plotting
    m_line = 0
    b_line = 0

    for i in range(ndim):
        mcmc = np.percentile(flat_samples[:, i], [16, 50, 84])
        q = np.diff(mcmc)
        #txt = "\mathrm{{{3}}} = {0:.3f}_{{-{1:.3f}}}^{{{2:.3f}}}"
        #txt = txt.format(mcmc[1], q[0], q[1], labels[i])
        print(labels[i], "=",mcmc[1], "+/-:",q[0],",", q[1])
        if i == 0: # CS m
            m_line = mcmc[1]
            print("m 16, 50, and 84th precentiles:",mcmc)
           
        if i == 1: # CS b
            b_line = mcmc[1]
            print("b 16, 50, and 84th precentiles:",mcmc)
            
    plt.show()
    
    plt.cla() # Clear plt so we don't overplot
    plt.scatter(x,y,color="black")
    plt.errorbar(x,y,yerr=yerr,color="black")
    x_lin = np.linspace(1,19)
    # CS Plot the model represented by our best-fitting parameters
    plt.scatter(x,m_line*x+b_line,color="blue")
    plt.plot(x_lin,m_line*x_lin+b_line,color="blue",linestyle="--")
    # CS Print best fit parameters to the plot
    bfit_lin = mpatches.Patch(facecolor='b', label=("m:"+str(round(m_line,3))+", b:" + str(round(b_line,3))), linewidth = 0.5, edgecolor = 'black')
    legend = plt.legend(handles=[bfit_lin], loc = 0, bbox_to_anchor = (1,1), fontsize = 10, fancybox = False, edgecolor = "none", facecolor = "none")
    plt.show()


def emcee_quad_fit(data_location):
    """ Uses a Bayseian framework to fit a quadratic
    model of the form y=a2*x**2+a1*x+a0 to the data 
    
    Parameters
    ----------
    data_location: :class: 'string file-path'
        File-path to data file with four columns, x,s,sig,yerr
    Returns
    -------
    A corner plot of parameters' posterior probability distributions
    A plot of y(x) including error bars
    Prints to screen best fitting parameters and uncertanties
    based on the 16th, 50th, and 84th precentiles
    """
    # CS Import dataxy
    data = Table.read(str(data_location))
    x = data["x"]
    y = data["y"]
    yerr = data["yerr"] # CS estimated error (sigma) on the y data
    def log_likelihood(theta, x, y, yerr):
        """ Computes the log likelihood of input parameters contained in theta
        for data parameterized by x, y and yerr
        
        Parameters
        ----------
        theta: :class: 'array'
            Input function parameters (i.e. m and b for a line)
        x: :class: 'array'
            Data x values
        y: :class: 'array'
            Data y values
        yerr: :class: 'array'
            Estimated error on the y data
        Returns
        -------
        Log liklihood for given parameters and data
        """
        a2, a1, a0  = theta
        model = a2 * x**2 + a1 * x + a0 # CS Quadratic model
        sigma2 = yerr**2 + model**2
        return -0.5 * np.sum((y - model) ** 2 / sigma2 + np.log(sigma2))
    
    # CS Choose some inital values that make sense for our line
    a2_true = -0.1
    a1_true = -0.5
    a0_true = 1

        
    np.random.seed(42)
    nll = lambda *args: -log_likelihood(*args)
    initial = np.array([a2_true,a1_true,a0_true]) + 0.1 * np.random.randn(3)
    soln = minimize(nll, initial, args=(x, y, yerr))
    a2_ml, a1_ml, a0_ml = soln.x
    
    
    # CS Make our prior bounds realistic for our line
    def log_prior(theta):
        """ Computes a uniform prior for an input array of parameters
        
        Parameters
        ----------
        theta: :class: 'array'
            Input function parameters (i.e. m and b for a line)
            
        Returns
        -------
        0 if parameters are within prior range, and -inf otherwise
        """
        a2, a1, a0  = theta
        if -5 < a2 < 10 and -6 < a1 < 3 and -1 < a0 < 10:
            # CS A uniform prior is zero when parameters are withinoru specified range
            return 0.0
        return -np.inf
    
    def log_probability(theta, x, y, yerr):
        """ Combines log_prior and log_likelihood functions
        to comput log probability
        
        Parameters
        ----------
        theta: :class: 'array'
            Input function parameters (i.e. m and b for a line)
        x: :class: 'array'
            Data x values
        y: :class: 'array'
            Data y values
        yerr: :class: 'array'
            Estimated error on the y data
        Returns
        -------
        log_probability as log_prior + log_likelihood
        """
        lp = log_prior(theta)
        if not np.isfinite(lp):
            return -np.inf
        return lp + log_likelihood(theta, x, y, yerr)
    
    # CS Start the walkers near the maximum likelihood
    pos = soln.x + 1e-4 * np.random.randn(32, 3)
    nwalkers, ndim = pos.shape
    
    sampler = emcee.EnsembleSampler(
        nwalkers, ndim, log_probability, args=(x, y, yerr)
    )
    # CS 5000 walkers should be fine according to emcee tutorials, check burn in though
    sampler.run_mcmc(pos, 5000, progress=True);
    
    # CS plot the time series of the parameters in the chain
    fig, axes = plt.subplots(3, figsize=(10, 7), sharex=True)
    samples = sampler.get_chain()
    labels = ["a2", "a1", "a0"]
    for i in range(ndim):
        ax = axes[i]
        ax.plot(samples[:, :, i], "k", alpha=0.3)
        ax.set_xlim(0, len(samples))
        ax.set_ylabel(labels[i])
        ax.yaxis.set_label_coords(-0.1, 0.5)
    
    axes[-1].set_xlabel("step number");
    
    # CS Print steps required to burn in
    tau = sampler.get_autocorr_time()
    print("Burn in:",tau)
    
    # CS Burn in takes a little less than 50 steps so 120 should be fine to discard
    # CS Remove the burn in steps and thin by 15 steps
    flat_samples = sampler.get_chain(discard=120, thin=15, flat=True)
    print("flat list of samples shape:",flat_samples.shape)
    # CS Note that removing the burn in steps lets the emcee walkers actually have time
    # CS to move to a "good" location in our parameter space before we start "trusting" the values
    # CS that they give us. 
    
    fig = corner.corner(flat_samples, labels=labels)
        
    # CS Get fit params for easy plotting
    a2 = 0
    a1 = 0
    a0 = 0

    for i in range(ndim):
        mcmc = np.percentile(flat_samples[:, i], [16, 50, 84])
        q = np.diff(mcmc)
        #txt = "\mathrm{{{3}}} = {0:.3f}_{{-{1:.3f}}}^{{{2:.3f}}}"
        #txt = txt.format(mcmc[1], q[0], q[1], labels[i])
        print(labels[i], "=",mcmc[1], "+/-:",q[0],",", q[1])
        # CS a2
        if i == 0:
            a2 = mcmc[1]
            print("a2 16, 50, and 84th precentiles:",mcmc)
        # CS a1
        if i == 1:
           a1 = mcmc[1]
           print("a1 16, 50, and 84th precentiles:",mcmc)
        # CS a0
        if i == 2:
            a0 = mcmc[1]
            print("a0 16, 50, and 84th precentiles:",mcmc)
           
    plt.show()
    
    plt.cla()
    plt.scatter(x,y,color="black")
    plt.errorbar(x,y,yerr=yerr,color="black")
    x_lin = np.linspace(1,19)
    # CS Plot the model represented by our best-fitting parameters
    plt.scatter(x,a2*x**2+a1*x+a0,color="green")
    plt.plot(x_lin,a2*x_lin**2+a1*x_lin+a0,color="green",linestyle="--")
    # CS Print best fit parameters to the plot
    bfit_quad = mpatches.Patch(facecolor='green', label=("a2:"+str(round(a2,3))+", a1:" + str(round(a1,3))+", a0:" +str(round(a0,3))), linewidth = 0.5, edgecolor = 'black')
    legend = plt.legend(handles=[bfit_quad], loc = 0, bbox_to_anchor = (1,1), fontsize = 10, fancybox = False, edgecolor = "none", facecolor = "none")
    plt.show()



if __name__ == "__main__":
    emcee_line_fit("/d/scratch/ASTR5160/final/dataxy.fits")
    emcee_quad_fit("/d/scratch/ASTR5160/final/dataxy.fits")
    print("-------")
    print("Comments:")
    print("-------")
    print("The posterior probability distribution for a2 shows a signifcinatly")
    print("narrower probability distribution (-0.15,0.15), with a defninitive spike and zero values")
    print("outside of the wings. On the contratry, in the line fit, m and b both show much")
    print("wider probability distributions (-1.6,0 for m), indicating a fit that is inferior to the quadratic fit")
    print("Thus we can say that the quadratic model is a better fit to our data than the linear model")
    print("(Saying: The quadratic model is the BEST fit... would be incorrect, we only know it is superior to a linear fit)")