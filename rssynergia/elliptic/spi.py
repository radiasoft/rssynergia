import sys
import synergia
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import gridspec


def get_coords(bunch):
    '''return basic coordinates'''
    
    xinits = np.abs(bunch[::,0]*1.e3) #initial x coordinates in mm by particle ID
    yinits = np.abs(bunch[::,2]*1.e3) #initial y coordinates in mm by particle ID
    xpinits = np.abs(bunch[::,1]*1.e3) #initial x' coordinates in mrad by particle ID
    ypinits = np.abs(bunch[::,3]*1.e3) #initial x' coordinates in mrad by particle ID
    
    return [xinits,yinits,xpinits,ypinits]


def gather_stats(hArray, numP):
    #Separate into H and I arrays
    h1 = hArray[:,:numP]
    h2 = hArray[:,numP:]
    #Note that h1, h2 are indexed by turn.
    #Thus, h1[turn] returns the array of h1 values across all macro particles for that turn, and h2[turn] does the same for I.
    hT1 = h1.transpose()
    hT2 = h2.transpose()
    #Note that hT1, hT2 are indexed by particle ID. 
    #Thus hT1[ID] returns the turn by turn array of H values for particle ID, and hT2[ID] does the same for I.

    #One approach is to plot the variation in H/I versus initial amplitude in the x-plane
    #Alternatively, we can plot mean H amplitude versus standard deviation in I.

    hmeans = []
    hstds = []
    imeans = []
    istds = []

    #We want to look at each particle individually, thus use hT arrays
    for partID in range(numP):
        h_mean = np.mean(hT1[partID]) * 1.e6
        h_std = 100*np.std(hT1[partID])*1.e6/h_mean
        i_mean = np.mean(hT2[partID]) * 1.e6
        i_std = 100*np.std(hT2[partID])*1.e6/i_mean
    
        #particle amplitudes in coordinate space
        #x_init = original_particles[partID,0]
        #xp_init = original_particles[partID,1]
        #y_init = original_particles[partID,2]
        #yp_init = original_particles[partID,3]
    
        hmeans.append(h_mean)
        imeans.append(i_mean)
        istds.append(i_std)
        hstds.append(h_std)
    
    
    return [hmeans,hstds,imeans,istds]
    
def print_stats(hArray, opts):
    
    '''Print a summary of cumulative H, I values to a text file in simulation output directory'''
    
    stats = gather_stats(hArray, opts.macro_particles)
    
    summary = [np.mean(stats[i]) for i in range(len(stats))]
    
    categories = ["H", "H_std", "I", "I_std"]
    
    
    
    if opts.printfile:
        #print stats to output file
        printfile = "{}/{}".format(opts.relpath,opts.printfile)
        
        with open(printfile, 'w') as f:
            f.write("Summary of invariants for {} particles with lattice {} \n".format(opts.macro_particles, opts.lattice_name))
            for index, val in enumerate(summary):
                if index in (0, 2):
                   f.write("{} [mm-mrad]: {} \n".format(categories[index], val))
                else:
                    f.write("{} [%]: {} \n".format(categories[index], val))

    #print to console regardless
    print"Summary of invariants for {} particles with lattice {}".format(opts.macro_particles, opts.lattice_name)
    for index, val in enumerate(summary):
        if index in (0, 2):
            print"{} [mm-mrad]: {}".format(categories[index], val)
        else:
            print"{} [%]: {}".format(categories[index], val)


def plot_H(coords, stats):
    
    '''Plot variations for H'''
    
    xinits = coords[0]
    yinits = coords[1]
    xpinits = coords[2]
    ypinits = coords[3]
    hstds = stats[1]
    
    from matplotlib import gridspec
    fig = plt.figure(figsize=(12,8))
    gs = gridspec.GridSpec(2, 2, width_ratios=[1,1], height_ratios = [1,1]) 

    ylabel = "$\sigma(H)$ %"

    ax0 = plt.subplot(gs[0])
    ax0.scatter(xinits,hstds)
    ax0.set_title("Initial $x$ versus $\sigma (H)$")
    ax0.set_xlabel("x [mm]")
    ax0.set_ylabel(ylabel)
    #ax0.set_aspect(aspect=2.0)

    ax1 = plt.subplot(gs[1])
    ax1.scatter(yinits,hstds)
    ax1.set_title("Initial $y$ versus $\sigma (H)$")
    ax1.set_xlabel("y [mm]")
    ax1.set_ylabel(ylabel)
    #ax1.set_aspect(aspect=2.0)

    ax2 = plt.subplot(gs[2])
    ax2.scatter(xpinits,hstds)
    ax2.set_title("Initial $x'$ versus $\sigma (H)$")
    ax2.set_xlabel("x' [mrad]")
    ax2.set_ylabel(ylabel)
    #ax2.set_aspect(aspect=2.0)

    ax3 = plt.subplot(gs[3])
    ax3.scatter(ypinits,hstds)
    ax3.set_title("Initial $y'$ versus $\sigma (H)$")
    ax3.set_xlabel("y' [mrad]")
    ax3.set_ylabel(ylabel)

    # Tweak spacing between subplots to prevent labels from overlapping
    #plt.subplots_adjust(hspace=2)

    #set figure title
    fig.canvas.set_window_title('Variation in H as a function of inital particle amplitude')
    fig.tight_layout()
    plt.show()    



def plot_I(coords, stats):

    '''Plot variations for I'''
    
    xinits = coords[0]
    yinits = coords[1]
    xpinits = coords[2]
    ypinits = coords[3]
    istds = stats[3]
    
    fig = plt.figure(figsize=(12,8))
    gs = gridspec.GridSpec(2, 2, width_ratios=[1,1], height_ratios = [1,1]) 

    ylabel = "$\sigma(I)$ %"

    ax0 = plt.subplot(gs[0])
    ax0.scatter(xinits,istds)
    ax0.set_title("Initial $x$ versus $\sigma (I)$")
    ax0.set_xlabel("x [mm]")
    ax0.set_ylabel(ylabel)
    #ax0.set_aspect(aspect=2.0)

    ax1 = plt.subplot(gs[1])
    ax1.scatter(yinits,istds)
    ax1.set_title("Initial $y$ versus $\sigma (I)$")
    ax1.set_xlabel("y [mm]")
    ax1.set_ylabel(ylabel)
    #ax1.set_aspect(aspect=2.0)

    ax2 = plt.subplot(gs[2])
    ax2.scatter(xpinits,istds)
    ax2.set_title("Initial $x'$ versus $\sigma (I)$")
    ax2.set_xlabel("x' [mrad]")
    ax2.set_ylabel(ylabel)
    #ax2.set_aspect(aspect=2.0)

    ax3 = plt.subplot(gs[3])
    ax3.scatter(ypinits,istds)
    ax3.set_title("Initial $y'$ versus $\sigma (I)$")
    ax3.set_xlabel("y' [mrad]")
    ax3.set_ylabel(ylabel)

    # Tweak spacing between subplots to prevent labels from overlapping
    #plt.subplots_adjust(hspace=2)

    #set figure title
    fig.canvas.set_window_title('Variation in I as a function of inital particle amplitude')
    fig.tight_layout()
    plt.show()


def plot_amplitudes(hArray, bunch, opts):
    
    coords = get_coords(bunch)
    stats = gather_stats(hArray,opts.macro_particles)
    num = opts.num
    
    
    if num == 1:
        plot_H(coords, stats)
    elif num == 2:
        plot_I(coords, stats)   
    else:
        print 'Invariant number specified as ' + str(num) + '. Index must be 1 or 2. Plotting 1st invariant stats.'
        plot_H(coords, stats)
