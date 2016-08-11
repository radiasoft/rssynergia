#Utilities for use with python data analysis, specifically for analyzing Synergia and Warp outputs.
#
#Author: Nathan Cook
# 10/28/2015

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt


def compare_distribution(arr1, arr2, numBins, norm=True):
    '''Construct and return a histogram comparing two distributions
    
    Arguments:
        arr1 (array-like): 1D array for binning
        arr2 (array-like): 1D array for binning
        numbins (int): # of bins for histogram
        
        (optional) norm: normalization flag (default True)
    
    Returns:
        fig (matplotlib.figure.Figure): - a matplotlib figure handle
    
    Constructs a matplotlib plot and returns it. Automatically displayed when using IPython backend.    
    
    '''

    myVals1, myBins1 = np.histogram(arr1,numBins)
    bincenters1 = 0.5*(myBins1[1:]+myBins1[:-1])
    bin_width1 = 1.0*(myBins1[1]-myBins1[0])
    #normalize
    myVals_norm1 = myVals1/(np.max(myVals1)*1.0)
    
    myVals2, myBins2 = np.histogram(arr2,numBins)
    bincenters2 = 0.5*(myBins2[1:]+myBins2[:-1])
    bin_width2 = 1.0*(myBins2[1]-myBins2[0])
    #normalize
    myVals_norm2 = myVals2/(np.max(myVals2)*1.0)    
    
    
    #Set some matplotlib standards for plots 
    mpl.rcParams['figure.figsize'] = 8, 6
    
    fig = plt.figure()
    ax = fig.gca()
    
    ax.set_xlabel("Array values", fontsize=14)
    ax.set_ylabel("Relative population", fontsize=14)
    
    if norm:
        ax.plot(bincenters1,myVals_norm1, c='r', label='1')
        ax.plot(bincenters2,myVals_norm2,  c='b', label='2')
    else:
        ax.plot(bincenters1,myVals1, c='r', label='1')
        ax.plot(bincenters2,myVals2,  c='b', label='2')        
    
    plt.legend(loc='best')
    
    #close the first display call
    plt.close()
    
    return fig #return the figure handle

    
    

def plot_distribution(arr, numBins, norm=False):
    '''Construct a histrogram, then plot the resulting distribution.
    
    Arguments:
        arr (array-like): 1D array for binning
        numbins (int): # of bins for histogram
        
        (optional) norm: normalization flag (default False)
    
    Returns:
        fig (matplotlib.figure.Figure): - a matplotlib figure handle
    
    Constructs a matplotlib plot and returns it. Automatically displayed when using IPython backend.
    
    '''
    myVals, myBins = np.histogram(arr,numBins)
    bincenters = 0.5*(myBins[1:]+myBins[:-1])
    bin_width = 1.0*(myBins[1]-myBins[0])
    #normalize
    myVals_norm = myVals/(np.max(myVals)*1.0)
    
    #Set some matplotlib standards for plots 
    mpl.rcParams['figure.figsize'] = 8, 6
    
    fig = plt.figure()
    ax = fig.gca()
    
    ax.set_xlabel("Array values", fontsize=14)
    ax.set_ylabel("Relative population", fontsize=14)
    
    ax.plot(bincenters,myVals_norm, c='k')
    #close the first display call
    plt.close()
    
    return fig #return the figure handle
    
    
def get_distribution(arr, numBins, norm = True):
    '''Construct a histrogram, and return the resulting distribution.
    
    Arguments:
        arr (array-like): 1D array for binning
        numbins (int): # of bins for histogram
        
        (optional) norm: normalization flag (default False)
    
    Returns:
        bincenters (ndarray): array of bincenter positions
        myVals (ndarray): array of histrogram values
    
    '''
    
    myVals, myBins = np.histogram(arr,numBins)
    bincenters = 0.5*(myBins[1:]+myBins[:-1])
    bin_width = 1.0*(myBins[1]-myBins[0])
    
    if norm:
        #normalize
        myVals_norm = myVals/(np.max(myVals)*1.0)
        return bincenters, myVals_norm
    else:
        #don't normalize
        return bincenters, myVals
        
import matplotlib as mpl

def plot_both(z, f_goal, f_actual):
    '''
    Plot two distributions with the same horizontal coordinates. Constructs a matplotlib 
    plot and returns it. Automatically displayed when using IPython backend.
        
    Arguments:
        z (ndarray): 1darray of floats (micron) - Array of the longitudinal coordinate shared between both dependant variables

        f_goal (ndarray): 1D array corresponding to values of the desired distribution at each z coordinate
        f_actual (ndarray): 1D array corresponding to values of the actual distribution at each z coordinate
    
    Returns:
        fig (matplotlib.figure.Figure): - a matplotlib figure handle
    
    Constructs a matplotlib plot and returns it. Automatically displayed when using IPython backend.
        
    '''
    
    #Set some matplotlib standards for plots 
    mpl.rcParams['figure.figsize'] = 8, 6
    
    fig = plt.figure()
    ax = fig.gca()
    
    ax.plot(z, f_goal, c='b', label = 'Goal Distribution')
    ax.plot(z, f_actual, c='g', label = 'Extracted Distribution')
    
    if min(z) < 0:
        min_z = 1.1*min(z)
    else:
        min_z = 0.9*min(z)
    
    if max(z) < 0:
        max_z = 0.9*max(z)
    else:
        max_z = 1.1*max(z) 
     
    ax.set_xlim([min_z,max_z])
    
    ax.set_xlabel("Longitudinal coordinate", fontsize=16)
    ax.set_ylabel("Relative population", fontsize=16)
    
    ax.set_title("Comparison of distributions", fontsize=16)
    
    plt.legend(loc='best')
    
    #close the display call
    plt.close()
    
    return fig #return the figure handle