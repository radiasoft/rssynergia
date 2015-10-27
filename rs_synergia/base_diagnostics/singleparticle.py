import sys
import os
import tables
import numpy as np
import math
import matplotlib
import matplotlib.pyplot as plt
import options
#from mpl_toolkits.axes_grid import make_axes_locatable
#from mpl_toolkits import axes_grid

coords = {}
coords['x'] = 0
coords['xp'] = 1
coords['y'] = 2
coords['yp'] = 3
coords['cdt'] = 4
coords['dpop'] = 5
coords['id'] = 6

#opts.coords = coords

fixed = {}
fixed['alpha_e'] = 1
fixed['beta_e'] = 1
fixed['beta'] = 0.6539
fixed['gamma'] = (1 + fixed['alpha_e']**2)/(fixed['beta_e'])
fixed['t'] = 0.4
fixed['c'] = 0.01

###################################### PLOTTERS ###################################

def plot_SPH(hamArray,opts, ID=0):
    
    '''
    Create plot for single particle hamiltonian for one particle over pre-specified turn #s
    
    Arguments:
    ham - numpy array of single particle hamiltonian values
    opts - options object
    
    Optional:
    ID - specify the specific particle to be plotted (default 0)
    
    '''
    
    #general plot settings
    matplotlib.rcParams['figure.autolayout'] = True
    
    #change ID only if necessary
    if opts.ID:
        ID = opts.ID
   
    h = np.arange(opts.turns+1) #plus 1 to account for initial conditions
    v = hamArray[:,ID]
    vinit = v[0] #take the iniital values as the normalization value
    
    if opts.norm:
        vScale = v/vinit
        ymin = 0
        ymax = 2
    else:
        vScale = v
        ymax = 1.5*vScale.max()
        ymin = 0.5*vScale.min()
    
    fig = plt.figure()
    plt.subplot(1,1,1)
    ax = plt.gca()
    
    #print ymin, ymax
    
    ax.scatter(h,vScale, c ='b', s=6)
    #ax.set_aspect('equal')
    
    ax.set_ylim([ymin,ymax])
    #plt.plot(h,v, 'o')
    plt.xlabel(opts.hcoord,fontsize=12)
    plt.ylabel(opts.vcoord,fontsize=12)
    title = 'Single particle hamiltonian for particle ' + str(ID)
    if not opts.lattice_name== None:
        title = title + ' for lattice ' + opts.lattice_name
    plt.title(title, y=1.05, fontsize=14)
    #plt.draw()
    #fig.tight_layout()
    plt.show()
    
    if opts.save:
        sv_title = 'SPH_0_'+ opts.lattice_name + '.pdf'
        fig.savefig(sv_title, bbox_inches='tight')


def plot_P(PArray, opts, num=10, ID=0):
    
    '''Create a Poincare plot for specified particle IDs over pre-specified turn #s
    
    
    Note: Currently, opts.plots must be a list of coordinates of length 2!
    
    '''

    #general plot settings
    matplotlib.rcParams['figure.autolayout'] = True
    
    if opts.num:
        num = opts.num
    
    #plot up to 10 particle tracks
    if PArray.shape[1] < num:
        num = PArray.shape[1]

    #plot specified # of turns
    if opts.turns:
        turns = opts.turns
    else: turns = PArray.shape[0]        
    
    plots = {'x':0, 'px':1, 'y':2, 'py':3}
    
    cropped = PArray[:turns,:num,(plots[opts.plots[0]],plots[opts.plots[1]])]
    
    reshaped = cropped.reshape(cropped.size/2,2).T
    
    #separate horizontal and vertical components
    h = reshaped[0]
    v = reshaped[1]
    
    #ymax = max(abs([h,v]))
    
    
    if opts.scale:
        fig = plt.figure(figsize=(opts.scale*8,opts.scale*6))
    else:
        fig = plt.figure(figsize=(8,6))
        
    plt.subplot(1,1,1)
    ax = plt.gca()
    
    #print ymin, ymax
    
    ax.scatter(h,v, c ='b', s=2)
    ax.set_aspect('equal') #want equal aspect ratios for Poincare plots
    
    #ax.set_ylim([ymin,ymax])
    #plt.plot(h,v, 'o')
    plt.xlabel(opts.hcoord,fontsize=round(12*opts.scale))
    plt.ylabel(opts.vcoord,fontsize=round(12*opts.scale))
    title = opts.hcoord + '-'+ opts.vcoord + ' for ' + str(turns) + ' turns'
    if not opts.lattice_name== None:
        title = title + ' for lattice ' + opts.lattice_name
    plt.title(title, y=1+0.05/opts.scale, fontsize=round(14*opts.scale))
    #plt.draw()
    #fig.tight_layout()
    plt.show()
    
    if opts.save:
        sv_title = 'Poincare'+'_' + opts.hcoord+'_' + opts.vcoord+'_'+ str(turns) + '_turns_'+  opts.lattice_name + '.pdf'
        fig.savefig(sv_title, bbox_inches='tight')     
    

def plot_J(JArray,opts, ID=0):
    
    '''
    Create plot for single particle invariant for one particle over pre-specified turn #s
    
    Arguments:
    JArray - numpy array of single particle hamiltonian values
    opts - options object
    
    Optional:
    ID - specify the specific particle to be plotted (default 0)
    
    '''
    
    #general plot settings
    matplotlib.rcParams['figure.autolayout'] = True
    
    #change ID only if necessary
    if opts.ID:
        ID = opts.ID
   
    if opts.turns+1 > JArray.shape[0]:
        print 'Not enough data to plot for {} turns. Plotting {} turns instead'.format(opts.turns,JArray.shape[0]-1)
        h = np.arange(JArray.shape[0])*opts.turnsPerDiag #plus 1 to account for initial conditions
        v = JArray[::,ID]
    else:
        h = np.arange(opts.turns+1)*opts.turnsPerDiag #plus 1 to account for initial conditions
        v = JArray[:opts.turns+1,ID]
    
    vinit = v[0] #take the iniital values as the normalization value
    
    if opts.norm:
        vScale = v/vinit
        ymin = 0
        ymax = 2
    else:
        vScale = v*1.e6 #convert to mm-mrad
        if opts.variance:
            ymax = (1+opts.variance)*vScale.mean()
            ymin = (1-opts.variance)*vScale.mean()
        else:
            ymax = 1.05*vScale.max()
            ymin = 0.95*vScale.min()            
    
    fig = plt.figure(figsize=(8,6))
    plt.subplot(1,1,1)
    ax = plt.gca()
    
    #print ymin, ymax
    
    ax.scatter(h,vScale, c ='b', s=6)
    #ax.set_aspect('equal')
    
    ax.set_ylim([ymin,ymax])
    #plt.plot(h,v, 'o')
    plt.xlabel(opts.hcoord,fontsize=12)
    plt.ylabel(opts.vcoord + " [mm-mrad]",fontsize=12)
    
    
    #title stuff
    if opts.num == 2:
        title = 'Second Invariant for particle ' + str(ID)
        sv_title = 'I_'+str(ID)+'_'+ opts.lattice_name + '.pdf'
    else:
        title = 'First invariant for particle ' + str(ID)
        sv_title = 'H_'+str(ID)+'_'+ opts.lattice_name + '.pdf'
        
    if not opts.elliptic:
        title = 'Courant synder invariant for particle ' + str(ID)
        sv_title = 'J_'+str(ID)+'_'+ opts.lattice_name + '.pdf'
        
    if not opts.lattice_name== None:
        title = title + ' for lattice ' + opts.lattice_name
        
    plt.title(title, y=1.05, fontsize=14)
    #plt.draw()
    #fig.tight_layout()
    plt.show()
    
    if opts.save:
        #sv_title = 'J_'+str(ID)+'_'+ opts.lattice_name + '.pdf'
        fig.savefig(sv_title, bbox_inches='tight') 


def plot_Both(array1,array2,opts, ID=0):
    
    '''plot 2 invariants on the same axis'''

    #general plot settings
    matplotlib.rcParams['figure.autolayout'] = True
    
    #change ID only if necessary
    if opts.ID:
        ID = opts.ID
   
    h = np.arange(opts.turns+1) #plus 1 to account for initial conditions
    v1 = array1[:opts.turns+1,ID]
    v2 = array2[:opts.turns+1,ID]
    vinit1 = v1[0] #take the iniital values as the normalization value
    vinit2 = v2[0] #
    
    if opts.norm:
        vScale1 = v1/vinit1
        vScale2 = v2/vinit2
        ymin = 0
        ymax = 2
    else:
        vScale1 = v1*1.e6 #convert to mm-mrad
        vScale2 = v2*1.e6
        if opts.variance:
            ymax = (1+opts.variance)*vScale1.mean()
            ymin = (1-opts.variance)*vScale1.mean()
        else:
            ymax = 1.05*vScale1.max()
            ymin = 0.95*vScale1.min()            
    
    if opts.scale:
        fig = plt.figure(figsize=(2*8,2*6))
    else:
        fig = plt.figure(figsize=(8,6))
    
    plt.subplot(1,1,1)
    ax = plt.gca()
    
    #print ymin, ymax
    
    ax.scatter(h,vScale1, c ='g', s=6)
    ax.scatter(h,vScale2, c = 'r', s=6)
    #ax.set_aspect('equal')
    
    ax.set_ylim([ymin,ymax])
    #plt.plot(h,v, 'o')
    plt.xlabel(opts.hcoord,fontsize=12)
    plt.ylabel(opts.vcoord + " [mm-mrad]",fontsize=12)
    
    
    #title stuff
    if opts.num == 2:
        title = 'Second Invariant for particle ' + str(ID)
        sv_title = 'I_'+str(ID)+'_'+ opts.lattice_name + '.pdf'
    else:
        title = 'First invariant for particle ' + str(ID)
        sv_title = 'H_'+str(ID)+'_'+ opts.lattice_name + '.pdf'
        
    if not opts.elliptic:
        title = 'Courant synder invariant for particle ' + str(ID)
        sv_title = 'J_'+str(ID)+'_'+ opts.lattice_name + '.pdf'
        
    if not opts.lattice_name== None:
        title = title + ' for lattice ' + opts.lattice_name
        
    plt.title(title, y=1.05, fontsize=14)
    #plt.draw()
    #fig.tight_layout()
    plt.show()
    
    if opts.save:
        #sv_title = 'J_'+str(ID)+'_'+ opts.lattice_name + '.pdf'
        fig.savefig(sv_title, bbox_inches='tight')     
    
    
#define an option to replicate the pltbunch.plot_bunch function?

def plot_tracks(tracks, opts, ID=0):
    '''
    
    Plot particle coordinates vs. 's' coordinate for a single particle.
    
    Arguments:
    tracks - array of coordinate values - organized to correspond to index of opts.plots!
    opts - options object
    
    Optional:
    ID - use to specify a particle - default None
    
    '''

    #general plot settings
    matplotlib.rcParams['figure.autolayout'] = True
    
    #change ID only if necessary
    if opts.ID:
        ID = opts.ID

    #if #turns specified, then slice on that as well
    if opts.turns:
        hcoord = np.arange(opts.turns+1)*opts.turnsPerDiag #plus 1 to account for initial conditions
        vcoords = tracks[:opts.turns+1,:,ID].transpose() 
    else:
        hcoord = np.arange(tracks.shape[0])*opts.turnsPerDiag 
        vcoords = tracks[::,:,ID].transpose()
   
    fig = plt.figure()
    plt.subplot(1,1,1)
    ax = plt.gca()
    
    mins = []
    maxes = []
    
    for index,v in enumerate(vcoords):
        #vinit = v[0] #take the iniital values as the normalization value
        scale = v.max() #take the max value to scale v
        
        if opts.norm:
            vScale = v/scale
            mins.append(0)
            maxes.append(2)
        else:
            vScale = v
            min1 = vScale.min()
            if min1 < 0:
                mins.append(1.5*min1)
            else:
                mins.append(0.5*min1)
            maxes.append(1.5*vScale.max())
        #ax.scatter(hcoord,vScale, c ='b', s=6)
        ax.plot(hcoord,vScale, label=opts.plots[index])
        
    ax.set_ylim([min(mins),max(maxes)])
    #plt.plot(h,v, 'o')
    plt.xlabel(opts.hcoord,fontsize=12)
    
    ylabel = ', '.join(opts.plots)
    plt.ylabel(ylabel,fontsize=12)
    title = 'Particle tracks for particle ' + str(ID)
    if not opts.lattice_name== None:
        title = title + ' for lattice ' + opts.lattice_name
    plt.title(title, y=1.05, fontsize=14)
    #plt.draw()
    #fig.tight_layout()
    plt.legend(loc='best')
    plt.show()
    
    if opts.save:
        sv_title = 'Tracks_' + str(ID) + '_' + opts.lattice_name + '.pdf'
        fig.savefig(sv_title, bbox_inches='tight')
    



################################## GETTERS ####################################

def get_particles(inputfile, lost=None):
    
    '''Reads an input file and returns a numpy array of particles and a dictionary of root values. 
    If lost particles specified, then returns those separately.
    
    '''
    
    f = tables.openFile(inputfile, 'r')
    particles = f.root.particles.read()
    
    #As a test, arbitrarily remove some particles - working as intended
    #particles = np.delete(particles,10,0)
    #lost.append(10)
    #particles = np.delete(particles,20,0)
    #lost.append(21)
    
    #define lost particles array
    #lost_particles = np.zeros((len(lost),7))
    #lost_particles[:,6] = lost
    
    #get appropriate reference properties from file root
    npart = particles.shape[0]
    mass = f.root.mass[()]
    p_ref = f.root.pz[()]
    sn = f.root.s_n[()] #period length
    tn = f.root.tlen[()] #cumulative tracked length for this file

    f.close()
    
    header = dict()
    header['n_part'] = npart
    header['mass'] = mass
    header['p_ref'] = p_ref
    header['s_val'] = sn
    header['t_len'] = tn
    
    #need to potentially adjust the counter!
    #adjustment = particles[npart-1, 6] - (npart-1)
    
    if lost:
    
        #define new lists
        lost_particles = []
        kept_particles = []
        lost_counter = list(lost)
    
        #separate lost particles
        for index,particle in enumerate(particles):
            #val = particle[6]-adjustment
            val = particle[6]
            if val in lost:
                lost_particles.append(particle)
                #remove from counter
                #lost_counter.remove(val)
            else:
                kept_particles.append(particle)
            
        #now we just need to make sure we fill out the kept_particles array to make it the proper length.
        #if not len(lost_particles) == len(lost):
        #    for num in lost_counter:
        #        #placeholder = [-1,-1,-1,-1,-1,-1,num]
        #        placeholder = [0, 0, 0, 0, 0, 0, num]
        #        lost_particles.append(placeholder)
    
        return header, np.asarray(kept_particles), np.asarray(lost_particles)
    
    else:
        #lost_particles = []
        
        return header, particles
    

def get_file_list(opts):

    '''
    
    Returns a list of files of the form 'particles*.h5' from the current directory 
    
    
    Optional Arguments:
    path - link to directory storing .h5 files - default None (executes within present working directory)
    
    '''
    
    if not opts.relpath:
        #If no relative path specified, check current directory 
        files = os.listdir(os.getcwd())
    else:
        #Otherwise check specified relative path
        path = os.path.join(os.getcwd(),opts.relpath)
        files = [os.path.join(path,fn) for fn in os.listdir(path)]
        #files = os.listdir(opts.path) #only change directories if a different path is specified

    pfiles = []

    #basic filtering for particles*.h5 files
    for filename in files:
            if filename.find('particles') > -1 and filename.endswith('.h5'):
                pfiles.append(filename)
    
    return pfiles
    
    
def get_lost_particle_list(opts):
    '''Returns a list of particle IDs corresponding to lost particles from inspection of the output files'''
    files = get_file_list(opts)
    
    #compare first file output to final file output
    header1, particles1 = get_particles(files[0])
    
    header2, particles2 = get_particles(files[-1])
    
    lost = []
    
    #first check if size is equal
    if not (header1['n_part'] == header2['n_part']):
        #make a boolean comparison
        boolVals = (particles1[0:header2['n_part'],6] == particles2[:,6])
        
        #construct list of lost particles
        boolList = list(boolVals)
        start = 0;
        stop = len(boolList);
        offset = 0; #offset must increment for each lost particle as list index dephases with particle index

        while start < stop:
            try:
                idx = boolList.index(False) +offset #include offset
                lost.append(int(particles1[idx,6])) #append particle ID, could be different from array index
                boolList.remove(False)
                start = idx
                offset +=1
                #print str(start) + 'is start and ' + str(stop) + 'is stop'
            except ValueError:
                start = stop
    
    #if only one particle is lost, then we have an int and not a list, so we want to cast as a list
    if type(lost) == int:
        return [lost]
    else: 
        return lost
    
    
def get_twiss(lattice_simulator):
    '''
    Returns an array of twiss parameters versus longitudinal coordinate 's' for a given lattice.
    
    Arguments:
    lattice_simulator - a Synergia lattice simulator
    
    Return values have array configuration: [s,betax,alphax,gammax,betay,alphay,gammay]
    '''
    
    lattice = lattice_simulator.get_lattice()
    
    twiss = []
    
    for elem in lattice.get_elements():
        temp = []
        
        lattice_functions=lattice_simulator.get_lattice_functions(elem)
        gamma_x = (1 + lattice_functions.alpha_x**2)/lattice_functions.beta_x
        gamma_y = (1 + lattice_functions.alpha_y**2)/lattice_functions.beta_y
        
        temp = [lattice_functions.arc_length, lattice_functions.beta_x, lattice_functions.alpha_x, gamma_x, lattice_functions.beta_y, lattice_functions.alpha_y, gamma_y]
        twiss.append(temp)
    
    #at the end, pre-pend final twiss values @ beginning with s-value 0
    twiss_init = twiss[-1]
    twiss_init[0]=0.0
    twiss.insert(0,twiss_init)
    
    #verify that we haven't inserted an improper 0 value at the end
    if twiss[-1][0] < twiss[-2][0]:
            twiss = twiss[:-1][:]
        
    return np.asarray(twiss)
    
    
def get_sliced_twiss(lattice_simulator):
    '''
    Returns an array of twiss parameters versus longitudinal coordinate 's' for a given lattice simulator.
    
    This represents an improvement in resolution of the lattice, as latticefunctions are calculated slice by slice. See lfplot.get_sliced_lattice_functions()
    
    Arguments:
    lattice_simulator - a Synergia lattice simulator
    
    Return values have array configuration: [s,betax,alphax,gammax,betay,alphay,gammay]
    '''
    
    lattice = lattice_simulator.get_lattice()
    
    twiss = []
    
    for aslice in lattice_simulator.get_slices():
        temp = []
        
        lattice_functions=lattice_simulator.get_lattice_functions(aslice)
        gamma_x = (1 + lattice_functions.alpha_x**2)/lattice_functions.beta_x
        gamma_y = (1 + lattice_functions.alpha_y**2)/lattice_functions.beta_y
        
        temp = [lattice_functions.arc_length, lattice_functions.beta_x, lattice_functions.alpha_x, gamma_x, lattice_functions.beta_y, lattice_functions.alpha_y, gamma_y]
        twiss.append(temp)
    
    #at the end, pre-pend final twiss values @ beginning with s-value 0
    twiss_init = twiss[-1]
    twiss_init[0]=0.0
    twiss.insert(0,twiss_init)
    
    #verify that we haven't inserted an improper 0 value at the end
    if twiss[-1][0] < twiss[-2][0]:
            twiss = twiss[:-1][:]
        
    return np.asarray(twiss)
    

def normalized_coordinates(header, particles, twiss, offset=None, units=None, ID=None):
    '''Return the an array of particles (fixed s) in normalized transverse coordinates rather than trace space
    
    Input Coordinates - x, x', y, y'
    Output Coordinates - sqrt(betax)*x, 
                        (1/sqrt(betax))*(betax*x' + alpha*x), 
                        sqrt(betay)*y, 
                        (1/sqrt(betay))*(betay*y + alpha*y)
    
    '''
    #consider offset for a few specific instances
    #if offset:
    #    sval = header['s_val'] + offset 
    #else:    
    #    sval = header['s_val']
    #sval = header['s_val'] + 4.9115191429
    #sval = header['s_val'] + 4.9115191429
    
    sval = header['s_val']
    
    #Must be careful to use the appropriate s_val, or else normalizaztion will be off
    svals = twiss[::,0]
    if sval > np.max(svals):
        sval = 0    #simply fix the s value to be 0 if the particle is past the final slice s-value.
    
    if offset:
        sval = sval+offset
    
    #print sval
    betaxvals = twiss[::,1]
    alphaxvals = twiss[::,2]
    gammaxvals = twiss[::,3]
    betayvals = twiss[::,4]
    alphayvals = twiss[::,5]
    gammayvals = twiss[::,6]
    
    #interpolate if needed
    if not sval in svals:
        betax = np.interp(sval, svals, betaxvals)
        betay = np.interp(sval, svals, betayvals)
        alphax = np.interp(sval, svals, alphaxvals)
        alphay = np.interp(sval, svals, alphayvals)
        gammax = np.interp(sval, svals, gammaxvals)
        gammay = np.interp(sval, svals, gammayvals)
    else:
        ind = list(svals).index(sval)
        betax = twiss[ind,1]
        alphax = twiss[ind,2]
        gammax = twiss[ind,3]
        betay = twiss[ind,4]
        alphay = twiss[ind,5]
        gammay = twiss[ind,6]
    
    #fix these values for center of NL element
    #betax = fixed['beta']
    #betay = betax
    #alphax = 0.0
    #alpahy = 0.0
    
    #betax = 0.5848554532
    #betay = betax
    #alphax = 0.0
    #alphay = 0.0
    
    
    #USED TO USE THIS
    #betax = fixed['beta']
    #betax = 6.1246858201346921
    #betay = betax
    #alphax = -1*3.0776835371752562
    #alphay = 1*alphax


    
    #alphax = -1*fixed['alpha_e']
    #alphay = -1*fixed['alpha_e']
    
    #if sval == 0:
    #    alphax = 1*fixed['alpha_e']
    #    alphay = 1*fixed['alpha_e']        
    
    #betax = 1
    #betay = 1
    
    #alphax = 0
    
    #print alphax
    #print alphay
    
    x = particles[:,coords['x']] #units m
    newx = x / math.sqrt(betax) #normalized
    
    xp = particles[:,coords['xp']] #unitless
    px = (betax*xp + alphax*x)/math.sqrt(betax) #normalized
    
    y = particles[:,coords['y']] #units m
    newy = y / math.sqrt(betay) #normalized
    
    yp = particles[:,coords['yp']] #unitless
    py = (betay*yp + alphay*y) /math.sqrt(betay) #normalized    
    
    #stack arrays then transpose to get array of coordinate vectors
    particles_norm = np.vstack((newx,px,newy,py)).T
    
    return particles_norm
    
def get_normalized_coords(filelist, twiss, lost=None, plotlost=False, num=None, ID=None):
    
    '''
    
    Returns a numpy array of normalized coordinate vectors obtained from .h5 files in filelist
    
    Arguments:
    filelist - A list of .h5 files containing particle array information
    twiss - array of twiss parameters for the lattice
    
    
    Returns a numpy array with dimensions #turns x #particles x #transverse coordinates(4).
    norms[B][A] returns vector of coordinates for particle A at turn B.
    '''

    #opts = options.Options()
    #opts.hcoord = 'x'
    #opts.vcoord = 'xp'
    #opts.lattice_name = 'FODO'
    norms = [] #norms is a list of arrays of particle coordinates
    
    for index,fileName in enumerate(filelist):
        inputfile = fileName
        
        if lost:
            header, particles, lost_particles = get_particles(inputfile, lost)
        else:
            header, particles = get_particles(inputfile)
        
        if plotlost:
            norm_coords = normalized_coordinates(header,lost_particles,twiss)
        else:
            norm_coords = normalized_coordinates(header, particles, twiss)
        
        #only append first num of tracked particles if not plotting all
        if num:
            norms.append(norm_coords[:num])
        else:  
            norms.append(norm_coords) #added flatten here
        #append 0s to maintain array size
        #if not header['n_part'] == norm_coords.shape[0]
    return np.asarray(norms)


def elliptic_coordinates(normalized, opts):
    '''Return the an array of elliptic coordinate values u and v for each 'particle' in the input array
    
    Input Coordinates - normal coordinates x, px, y, py
    Output Coordinates - u, v - vectors of length equal to the # of particles
    
    Arguments:
        - t, c - nonlinear magnet parameters (via opts)
    
    '''
    #beta = 0.6539
    #beta =  fixed['beta']
    beta = 6.1246858201346921
    
    t = opts.t 
    #c = opts.c 
    c = opts.c / np.sqrt(beta)
        
    x = normalized[:,0]
    #px = normalized[:,1]
    y = normalized[:,2]
    #py = normalized[:,3]
    #first need to adjust x and y by the c factor
    
    x = x/c
    y = y/c
    
    #this needs to be adjusted so that I work on the entire array in one swoop
    
    u = 0.5*(np.sqrt((x + 1.)**2 + y**2) + np.sqrt((x -1.)**2 + y**2))
    v = 0.5*(np.sqrt((x + 1.)**2 + y**2) - np.sqrt((x -1.)**2 + y**2))
    
    #f2u = u * np.sqrt(u**2 -1.) * np.arccosh(u)
    #g2v = v * np.sqrt(1. - v**2) * (-np.pi/2 + np.arccos(v))
    
    
    #elliptic = f2u + g2v / (u**2 - v**2)
    
    return [u,v]
    
    
def elliptic_hamiltonian(u,v, opts):
    
    '''
    
    Returns arrays of values for the first elliptic invariant (Hamiltonian) for a system with NL magnetic potential
    
    '''
    #beta = 0.6539
    beta =  fixed['beta']
    beta = 6.1246858201346921
    
    t = opts.t 
    #
    c = opts.c / np.sqrt(beta)
    #c = opts.c * np.sqrt(beta)
    
    f2u = u * np.sqrt(u**2 -1.) * np.arccosh(u)
    g2v = v * np.sqrt(1. - v**2) * (-np.pi/2 + np.arccos(v))
    
    elliptic = (f2u + g2v) / (u**2 - v**2)
    kfac = -1*t*c*c
    
    return kfac*elliptic
    
    


def second_invariant(normalized, u,v, opts):

    '''
    
    Returns an array containing the 2nd invariant for IOTA particles with the NL magnet
    
    Arguments:
    normalized - array of particles normalized coordinates
    u,v - arrays of elliptic coordinates
    c, t - elliptic potential strength parameters (via opts)
    
    '''
    
    #beta = 0.6539
    beta =  fixed['beta']
    beta = 6.1246858201346921
    
    t = -1*opts.t
    #c = opts.c
    #c = opts.c * np.sqrt(beta)
    c = opts.c / np.sqrt(beta)
    
    x = normalized[:,0]
    px = normalized[:,1]
    y = normalized[:,2]
    py = normalized[:,3]
    
    
    p_ang = (x*py - y*px)**2
    p_lin = (px*c)**2
    
    #harmonic part of potential
    f1u = c**2 * u**2 * (u**2 -1.)
    g1v = c**2 * v**2 * (1.-v**2)
    
    #elliptic part of potential
    f2u = -t * c**2 * u * np.sqrt(u**2-1.) * np.arccosh(u)
    g2v = -t * c**2 * v * np.sqrt(1.-v**2) * (0.5*np.pi - np.arccos(v))
    
    #combined - Adjusted this from Stephen's code
    fu = (0.5 * f1u - f2u)
    gv = (0.5 * g1v + g2v)
    
    #This is Stephen's previous statement, not sure if its correct
    #fu = 0.5 * f1u - f2u
    #gv = 0.5 * g1v + g2v
    
    invariant = (p_ang + p_lin) + 2.*(c**2) * (fu * v**2 + gv * u**2)/(u**2 - v**2)
    
    return invariant
       

def get_single_particle_elliptic_invariants(filelist, twiss, opts, lost, num=1):

    '''
    
    Returns an array of single particle invariants (Hamiltonian with elliptic potential)
    
    Arguments:
    filelist - list of output (.h5) files
    twiss - an array of twiss functions for the lattice
    opts - list of options (includes elliptic potential scalings c,t)
    num - 1 for 1st invariant, 2 for 2nd invariant

    
    '''

    #print num
    invariant = [] #invariant is a list of arrays of macroparticles
    
    for index,fileName in enumerate(filelist):
        inputfile = fileName
        
        if lost:
            header, particles, lost_particles = get_particles(inputfile, lost)
        else:
            header, particles = get_particles(inputfile, lost)
        
        #normalize coordinates    
        norm_coords = normalized_coordinates(header, particles, twiss)
        #construct elliptic coordinates for that file
        u,v = elliptic_coordinates(norm_coords, opts)
        
        
        if num == 1:
            #calculate the first invariant
            vals = single_particle_hamiltonian(norm_coords) + elliptic_hamiltonian(u,v,opts)            
        elif num == 2:
            #calculate the second invariant
            vals = second_invariant(norm_coords, u,v, opts)
        else:
            #otherwise calculate the first invariant
            print "Improper invariant number specified. Calculating 1st invariant by default."
            vals = single_particle_hamiltonian(norm_coords) + elliptic_hamiltonian(u,v,opts)
        
        invariant.append(vals)
        
        #invariant.append(single_particle_invariant(header, particles, twiss))

    return np.asarray(invariant)

def get_toy_twiss(vals, lat):
    
    '''Returns a toy twiss array containing fixed values for constructed 's' values '''
    
    #reiterate
    fixed['beta_e'] = vals[1]
    fixed['alpha_e'] = vals[2]
    fixed['beta'] = vals[1] #just make both of them that for now
    fixed['gamma'] = (1 + fixed['alpha_e']**2)/(fixed['beta_e'])
    
    svals = np.linspace(0,lat.get_length(),100)
    
    length = len(svals)
    twiss = np.zeros((length,7))    
    
    twiss[::,0] = svals
    twiss[::,1] = [fixed['beta'] for i in range(length)]
    twiss[::,2] = [fixed['alpha_e'] for i in range(length)]
    twiss[::,3] = [fixed['gamma'] for i in range(length)] 
    twiss[::,4] = [fixed['beta'] for i in range(length)]
    twiss[::,5] = [fixed['alpha_e'] for i in range(length)]
    twiss[::,6] = twiss[::,3]

    return twiss

def single_particle_invariant(header, particles, twiss, units=None, ID=None):
    
    '''
    
    Returns an array of single particle invariants (Courant Synder / Hamiltonian)
    
    Arguments:
    header - a header dictionary obtained from 'get_particles()'
    particles - an array of particles obtained from .h5 files via f.root.particles.read()
    twiss - an array of twiss functions
    
    Optional:
    ID - use to specify values for a single particle - default None
    
    '''
    
    sval = header['s_val']
    svals = twiss[::,0]
    betaxvals = twiss[::,1]
    alphaxvals = twiss[::,2]
    gammaxvals = twiss[::,3]
    betayvals = twiss[::,4]
    alphayvals = twiss[::,5]
    gammayvals = twiss[::,6]
    
    #interpolate if needed
    if not sval in svals:
        betax = np.interp(sval, svals, betaxvals)
        betay = np.interp(sval, svals, betayvals)
        alphax = np.interp(sval, svals, alphaxvals)
        alphay = np.interp(sval, svals, alphayvals)
        gammax = np.interp(sval, svals, gammaxvals)
        gammay = np.interp(sval, svals, gammayvals)
    else:
        ind = list(svals).index(sval)
        betax = twiss[ind,1]
        alphax = twiss[ind,2]
        gammax = twiss[ind,3]
        betay = twiss[ind,4]
        alphay = twiss[ind,5]
        gammay = twiss[ind,6]
    
    
    x = particles[:,coords['x']] #units m
    xp = particles[:,coords['xp']] #unitless
    
    y = particles[:,coords['y']] #units m
    yp = particles[:,coords['yp']] #unitless     
    
    invariantx = gammax*x**2 + 2*alphax*x*xp + betax*xp**2
    invarianty = gammay*y**2 + 2*alphay*y*yp + betay*yp**2
    
    #inv = invariantx**2+invarianty**2
    inv2 = invariantx + invarianty
    
    return inv2
    
def get_invariants(filelist, twiss, lost):
    
    '''
    
    Returns a numpy array of single particle invariants (Courant Synder) values obtained from .h5 files in filelist
    
    Arguments:
    filelist - A list of .h5 files containing particle array information
    twiss - array of twiss parameters for the lattice
    
    '''

    #opts = options.Options()
    #opts.hcoord = 'x'
    #opts.vcoord = 'xp'
    #opts.lattice_name = 'FODO'
    invariant = [] #hamiltonian is a list of arrays of macroparticles
    
    for index,fileName in enumerate(filelist):
        inputfile = fileName
        
        if lost:
            header, particles, lost_particles = get_particles(inputfile, lost)
        else:
            header, particles = get_particles(inputfile, lost)
            
        norm_coords = normalized_coordinates(header, particles, twiss)
        invariant.append(single_particle_hamiltonian(norm_coords))
        #invariant.append(single_particle_invariant(header, particles, twiss))
        
    return np.asarray(invariant)


def single_particle_hamiltonian(normalized, ID=None):
    
    '''
    
    UPDATED: Returns the single particle hamiltonian (quadratic component) in absence of nonlinearities 
    
    Arguments:
    normalized - the normalized coordinates for the particles
    
    Optional:
    ID - use to specify values for a single particle - default None
    
    '''
    
    #pref = header['p_ref'] #reference momentum in GeV/c
    #mass = header['mass'] #mass in GeV/c^2
    #gevc = 5.34428576e-19 #mkg/s per GeV/c
    
    x = normalized[:,0]
    px = normalized[:,1]
    y = normalized[:,2]
    py = normalized[:,3]
    
    #x = particles[:,coords['x']] #units m
    #px = particles[:,coords['xp']] #unitless
    
    #y = particles[:,coords['y']] #units m
    #py = particles[:,coords['yp']] #unitless 
    
    #if units:
    #    print "Units flag specified"
    #    px = px*pref*gevc #units kgm/s
    #    py = py*pref*gevc #units kgm/s
    
    if ID:
        #quadratic is the basis for calculating the hamiltonian, in absence of nonlinear couplings
        quadratic = 0.5* (px[ID]**2 + py[ID]**2) + 0.5*(x[ID]**2 + y[ID]**2)
    else:
        quadratic = 0.5* (px**2 + py**2) + 0.5*(x**2 + y**2)

    return quadratic
        
    

def get_hamiltonians(filelist):
    
    '''
    
    DEPRECATED!
    
    Returns a numpy array of single particle hamiltonian values obtained from .h5 files in filelist
    
    Arguments:
    filelist - A list of .h5 files containing particle array information
    
    '''

    #opts = options.Options()
    #opts.hcoord = 'x'
    #opts.vcoord = 'xp'
    #opts.lattice_name = 'FODO'
    hamiltonian = [] #hamiltonian is a list of arrays of macroparticles
    
    for index,fileName in enumerate(filelist):
        inputfile = fileName
        header, particles = get_particles(inputfile)
        hamiltonian.append(single_particle_hamiltonian(header, particles))
        
    return np.asarray(hamiltonian)

def get_tracks(filelist, opts):
    
    '''
    
    Returns a numpy array of single particle coordinate values obtained from .h5 files in filelist
    
    Arguments:
    filelist - A list of .h5 files containing particle array information
    
    '''

    #opts = options.Options()
    #opts.hcoord = 'x'
    #opts.vcoord = 'xp'
    #opts.lattice_name = 'FODO'
    tracks = [] #hamiltonian is a list of arrays of macroparticles
    
    for index,fileName in enumerate(filelist):
        inputfile = fileName
        header, particles = get_particles(inputfile)
        vec = []
        for coord in opts.plots:
            #make sure specified plot options are obtainable
            assert coord in coords.keys(), "Specified plot, %s is not available from: %s" %(coord, coords.keys())   
            vec.append(particles[:,coords[coord]])
        tracks.append(vec)
        
    return np.asarray(tracks)    


#######################################################################################################################################


def plot_Poincare(opts, noTwiss=False):
    
    '''Plot a poincare section in the desired normalized coordinates'''
    
    opts.hcoord = opts.plots[0]
    opts.vcoord = opts.plots[1]
    
    files = get_file_list(opts)
    lost = get_lost_particle_list(opts)
    #twiss = get_twiss(opts.lattice_simulator)
    #if noTwiss:
    #    twiss = noTwiss
    #else:
    twiss = get_sliced_twiss(opts.lattice_simulator)
    
    if opts.plot_lost:
         pArray = get_normalized_coords(files,twiss,lost,True)
        
    else:
        pArray = get_normalized_coords(files,twiss,lost)
    
    plot_P(pArray, opts) 
    
def toy_plot_Poincare(opts, vals):
    
    '''Plot a poincare section in the desired normalized coordinates for the toy R-matrix simulations'''
    
    opts.hcoord = opts.plots[0]
    opts.vcoord = opts.plots[1]
    
    files = get_file_list(opts)
    lost = get_lost_particle_list(opts)

    fixed['beta_e'] = vals[1]
    fixed['alpha_e'] = vals[2]
    fixed['beta'] = vals[1] #just make both of them that for now
    fixed['gamma'] = (1 + fixed['alpha_e']**2)/(fixed['beta_e'])

    twiss = get_toy_twiss(vals, opts.lattice)
    
    
    if opts.plot_lost:
         pArray = get_normalized_coords(files,twiss,lost,True)
        
    else:
        pArray = get_normalized_coords(files,twiss,lost)
    
    plot_P(pArray, opts) 


def plot_Invariant(opts):
    
    '''
    
    Plots the single particle hamiltonian over a specified # of turns and # of particles
    
    Arguments:
    opts - an Options object specifying # of turns, particle #s, etc.
    
    '''
    
    opts.hcoord = 'turn #'
    opts.vcoord = 'J(p,q)'
    opts.elliptic = False

    files = get_file_list(opts)
    #twiss = get_twiss(opts.lattice_simulator)
    twiss = get_sliced_twiss(opts.lattice_simulator)
    lost = get_lost_particle_list(opts)
    jArray = get_invariants(files, twiss, lost)
    #return jArray
    plot_J(jArray,opts)
    
def plot_elliptic_Invariant(opts):
    '''
    
    Plots the single particle hamiltonian for NL elliptic potential over a specified # of turns and # of particles
    
    Arguments:
    opts - an Options object specifying # of turns, particle #s, etc.
    
    '''
    
    opts.hcoord = 'turn #'
    #opts.t = 0.4
    #opts.c = 0.01
    opts.t = fixed['t']
    opts.c = fixed['c']
    opts.elliptic = True
    
    if opts.num:
        num = opts.num
    else:
        num = 1
        opts.num = 1
        
    if num == 1:
        opts.vcoord = 'H(p,q)'
    else:
        opts.vcoord = 'I(p,q)'
    
    files = get_file_list(opts)
    twiss = get_sliced_twiss(opts.lattice_simulator)
    #t2 = twiss[:-1,:]
    #twiss = get_twiss(opts.lattice_simulator)
    lost = get_lost_particle_list(opts)
    jArray = get_single_particle_elliptic_invariants(files, twiss, opts, lost, num)
    #jArray = get_invariants(files, twiss, lost)
    #return jArray
    plot_J(jArray,opts)
    
    
    
    

def toy_plot_elliptic_Invariant(opts, vals):
    '''
    
    Plots the single particle hamiltonian for NL elliptic potential over a specified # of turns and # of particles
    
    Adapted for toy model with R-matrix
    
    Arguments:
    opts - an Options object specifying # of turns, particle #s, etc.
    
    '''
    
    opts.hcoord = 'turn #'
    #opts.t = 0.4
    #opts.c = 0.01
    opts.t = fixed['t']
    opts.c = fixed['c']
    opts.elliptic = True
    
    if opts.num:
        num = opts.num
    else:
        num = 1
        opts.num = 1
        
    if num == 1:
        opts.vcoord = 'H(p,q)'
    else:
        opts.vcoord = 'I(p,q)'
    
    files = get_file_list(opts)
    
    fixed['beta_e'] = vals[1]
    fixed['alpha_e'] = vals[2]
    fixed['beta'] = vals[1] #just make both of them that for now
    fixed['gamma'] = (1 + fixed['alpha_e']**2)/(fixed['beta_e'])

    twiss = get_toy_twiss(vals, opts.lattice)
    #t2 = twiss[:-1,:]
    #twiss = get_twiss(opts.lattice_simulator)
    lost = get_lost_particle_list(opts)
    jArray = get_single_particle_elliptic_invariants(files, twiss, opts, lost, num)
    #jArray = get_invariants(files, twiss, lost)
    #return jArray
    plot_J(jArray,opts)
    

def calc_elliptic_Invariant(opts):
    '''
    
    Returns the single particle hamiltonian for NL elliptic potential over a specified # of turns and # of particles
    
    Arguments:
    opts - an Options object specifying # of turns, particle #s, etc.
    
    '''
    
    opts.hcoord = 'turn #'
    #opts.t = 0.4
    #opts.c = 0.01
    opts.t = fixed['t']
    opts.c = fixed['c']
    opts.elliptic = True
    
    if opts.num:
        num = opts.num
    else:
        num = 1
        opts.num = 1
        
    if num == 1:
        opts.vcoord = 'H(p,q)'
    else:
        opts.vcoord = 'I(p,q)'
    
    files = get_file_list(opts)
    twiss = get_sliced_twiss(opts.lattice_simulator)
    #t2 = twiss[:-1,:]
    #twiss = get_twiss(opts.lattice_simulator)
    lost = get_lost_particle_list(opts)
    hArray = get_single_particle_elliptic_invariants(files, twiss, opts, lost, 1)
    iArray = get_single_particle_elliptic_invariants(files, twiss, opts, lost, 2)
    #jArray = get_invariants(files, twiss, lost)
    #return jArray
    return np.hstack((hArray,iArray))
    
    

def toy_calc_elliptic_Invariant(opts, vals):
    '''
    
    Returns the single particle hamiltonian for NL elliptic potential over a specified # of turns and # of particles
    
    Arguments:
    opts - an Options object specifying # of turns, particle #s, etc.
    
    '''
    
    opts.hcoord = 'turn #'
    #opts.t = 0.4
    #opts.c = 0.01
    opts.t = fixed['t']
    opts.c = fixed['c']
    opts.elliptic = True
    
    if opts.num:
        num = opts.num
    else:
        num = 1
        opts.num = 1
        
    if num == 1:
        opts.vcoord = 'H(p,q)'
    else:
        opts.vcoord = 'I(p,q)'
        
        
    fixed['beta_e'] = vals[1]
    fixed['alpha_e'] = vals[2]
    fixed['beta'] = vals[1] #just make both of them that for now
    fixed['gamma'] = (1 + fixed['alpha_e']**2)/(fixed['beta_e'])
    
    files = get_file_list(opts)
    twiss = get_toy_twiss(vals, opts.lattice)
    #t2 = twiss[:-1,:]
    #twiss = get_twiss(opts.lattice_simulator)
    lost = get_lost_particle_list(opts)
    hArray = get_single_particle_elliptic_invariants(files, twiss, opts, lost, 1)
    iArray = get_single_particle_elliptic_invariants(files, twiss, opts, lost, 2)
    #jArray = get_invariants(files, twiss, lost)
    #return jArray
    return np.hstack((hArray,iArray))

def toy_stats_Invariant(hArray, opts):
    '''
    
    Return statistics on single particle invariants
    
    Arguments:
    hArray - horizontally stacked - 1st invariant is hArray[:,0], 2nd invariant is hArray[:,-1]
    
    '''
    array1 = hArray[:,opts.ID]
    array2 = hArray[:,opts.ID + opts.macro_particles]
    
    h_mean = np.mean(array1) * 1.e6
    h_std = np.std(array1) * 1.e6 
    i_mean = np.mean(array2) * 1.e6
    i_std = np.std(array2) * 1.e6
    
    print "H -  Mean: " + str(h_mean) + " [mm-mrad] std (%): " + str(100*h_std/h_mean)
    print "I -  Mean: " + str(i_mean) + " [mm-mrad] std (%): " + str(100*i_std/i_mean)


def stats_Invariant(hArray, opts):
    '''
    
    Return statistics on single particle invariants
    
    Arguments:
    hArray - horizontally stacked - 1st invariant is hArray[:,0], 2nd invariant is hArray[:,-1]
    
    '''
    array1 = hArray[:,opts.ID]
    array2 = hArray[:,opts.ID + opts.macro_particles]
    
    h_mean = np.mean(array1) * 1.e6
    h_std = np.std(array1) * 1.e6 
    i_mean = np.mean(array2) * 1.e6
    i_std = np.std(array2) * 1.e6
    
    print "H -  Mean: " + str(h_mean) + " [mm-mrad] std (%): " + str(100*h_std/h_mean)
    print "I -  Mean: " + str(i_mean) + " [mm-mrad] std (%): " + str(100*i_std/i_mean)
       

def plot_H_I(opts):
    '''
    
    Plots both the quadratic and elliptic components of the potential on the same axes. Useful for diagnosing coordinate transforms.
    
    '''
    
    opts.hcoord = 'turn #'
    #opts.t = 0.4
    #opts.c = 0.01
    opts.t = fixed['t']
    opts.c = fixed['c']
    opts.elliptic = True
    
        
    opts.vcoord = 'H(p,q), I(p,q)'
    
    files = get_file_list(opts)
    twiss = get_sliced_twiss(opts.lattice_simulator)
    lost = get_lost_particle_list(opts)
    
    num = 1
    IArray = get_single_particle_elliptic_invariants(files, twiss, opts, lost, num)
    
    HArray = get_invariants(files, twiss, lost)
    #jArray = get_invariants(files, twiss, lost)
    #return jArray
    #plot_J(jArray,opts)    
    plot_Both(HArray,IArray,opts)


def plot_Hamiltonian(opts):
    
    '''
    DEPRECATED
    
    Plots the single particle hamiltonian over a specified # of turns and # of particles
    
    Arguments:
    opts - an Options object specifying # of turns, particle #s, etc.
    
    '''
    
    opts.hcoord = 's'
    opts.vcoord = 'H(p,q)'
    
    files = get_file_list(opts)
    hamArray = get_hamiltonians(files)
    plot_SPH(hamArray,opts)
    
    
def track(opts):
    
    '''
    Plots specified coordinates for a single particle versus longitudinal coordinate `s`.
    
    Arguments:
    opts - an Options object specifying # of coordinates, turns, particle #s, etc.
    
    '''
    
    opts.hcoord = 'turn #'
    files = get_file_list(opts)
    tracks = get_tracks(files, opts)
    plot_tracks(tracks, opts)