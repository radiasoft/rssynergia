#Module for performing analysis on particle files diagnostics and bunch objects in Synergia.
#Used for plotting and statistical analysis of invariants and particle coordinates.
#Author: Nathan Cook
#4/21/2016


import sys
import os
import tables
import itertools
import numpy as np
import math
import matplotlib
import matplotlib.pyplot as plt
import options
import synergia

coords = {}
coords['x'] = 0
coords['xp'] = 1
coords['y'] = 2
coords['yp'] = 3
coords['cdt'] = 4
coords['dpop'] = 5
coords['id'] = 6

###################################### PLOTTERS ###################################


def plot_P(PArray, opts, num=10, ID=0):
    
    '''Create a Poincare plot for specified particle IDs over pre-specified turn #s.
    
    Arguments:
        PArray (ndArray): An array of particle data
        opts (options.Options): A Synergia options instance including an opts.plots list with >= 2 coordinates
        num (Optional[int]): The number of particles whose coordinates are plotted. Defaults to 10.
        ID (Optional[int]): The particle ID to plot when looking at a single particle. Defaults to 0.
    
    
    This function will save the plot if opts.save = True.
    
    '''

    #general plot settings
    matplotlib.rcParams['figure.autolayout'] = True
    
    if opts.num:
        num = opts.num

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
    
    if opts.scale:
        fig = plt.figure(figsize=(opts.scale*8,opts.scale*6))
    else:
        fig = plt.figure(figsize=(8,6))
        
    plt.subplot(1,1,1)
    ax = plt.gca()
    
    ax.scatter(h,v, c ='b', s=2)
    ax.set_aspect('equal') #want equal aspect ratios for Poincare plots

    plt.xlabel(opts.hcoord,fontsize=round(12*opts.scale))
    plt.ylabel(opts.vcoord,fontsize=round(12*opts.scale))
    title = opts.hcoord + '-'+ opts.vcoord + ' for ' + str(turns) + ' turns'
    if not opts.lattice_name== None:
        title = title + ' for lattice ' + opts.lattice_name
    if opts.plot_lost:
            title = 'Lost Particles: ' + title
    plt.title(title, y=1+0.05/opts.scale, fontsize=round(14*opts.scale))
    plt.show()
    
    if opts.save:
        sv_title = 'Poincare'+'_' + opts.hcoord+'_' + opts.vcoord+'_'+ str(turns) + '_turns_'+  opts.lattice_name + '.pdf'
        fig.savefig(sv_title, bbox_inches='tight')     
    

def plot_J(JArray,opts, ID=0):
    
    '''
    Create plot for single particle invariant for one particle over pre-specified turn numbers.
    
    Arguments:
        JArray (ndArray): An array of particle invariant data.
        opts (options.Options): A Synergia options instance including an opts.turns value.
        ID (Optional[int]): The particle ID to plot when looking at a single particle. Defaults to 0. 
                            Can be overwritten by specifying opts.ID
    
    This function will save the plot if opts.save = True.
    
    '''
    
    #general plot settings
    matplotlib.rcParams['figure.autolayout'] = True
    
    #change ID only if necessary
    if opts.ID:
        ID = opts.ID
   
    if opts.turns+1 > JArray.shape[0]:
        print 'Not enough data to plot for {} turns. Plotting {} turns instead'.format(opts.turns,JArray.shape[0]-1)
        h = np.arange(JArray.shape[0]) #plus 1 to account for initial conditions
        v = JArray[::,ID]
    else:
        h = np.arange(opts.turns+1) #plus 1 to account for initial conditions
        v = JArray[:opts.turns+1,ID]
    
    vinit = v[0] #take the iniital values as the normalization value
    
    if opts.norm:
        vScale = v/vinit
        ymin = 0
        ymax = 2
    else:
        if opts.num == 2:
            vScale = v*1.e9 #convert to 10^-9 m^2-rad^2 for second invariant
        else:
            vScale = v*1.e6 #convert to mm-mrad for all others
        if opts.variance:
            ymax = (1+opts.variance)*vScale.mean()
            ymin = (1-opts.variance)*vScale.mean()
        else:
            ymax = 1.05*vScale.max()
            ymin = 0.95*vScale.min()            
    
    fig = plt.figure(figsize=(8,6))
    plt.subplot(1,1,1)
    ax = plt.gca()

    ax.scatter(h,vScale, c ='b', s=6)

    ax.set_ylim([ymin,ymax])
    plt.xlabel(opts.hcoord,fontsize=12)
    
    if opts.num == 2:
        plt.ylabel(opts.vcoord + " [10$^{-9}$ m$^2$-rad$^2$]",fontsize=12)
    else:
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

    plt.show()
    
    if opts.save:
        fig.savefig(sv_title, bbox_inches='tight') 



################################## GETTERS ####################################

def get_one_particle(inputfile, ID=1): 
    '''
    Reads an input file and returns a single particle's coordinates specified by particle ID.
    
    Arguments:
        inputfile (str): path to a .h5 file containing particle diagnostics.
        ID (Optional[int]): particle ID for the one particle to get. Defaults to 1.
        
    Returns:
        particle (ndArray): array of particle data [x, x', y, y', cdt, dp, ID]
    
    '''
    
    f = tables.openFile(inputfile, 'r')
    particles = f.root.particles.read()
    
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
    
    #separate lost particles
    for particle in particles:
        val = particle[6]
        if val == ID:
            return particle

def get_some_particles(inputfile, ID=np.arange(100)):
    '''
    Reads an input file and returns a several particles' coordinates specified by particle ID.
    
    Arguments:
        inputfile (str): path to a .h5 file containing particle diagnostics.
        ID (Optional[list]): list of particle ID for the one particle to get. Defaults to [1:100]
        
    Returns:
        part_vals (ndArray): array of particle data [x, x', y, y', cdt, dp, ID] for each particle
    
    '''
    
    f = tables.openFile(inputfile, 'r')
    particles = f.root.particles.read()
    
    #get appropriate reference properties from file root
    npart = particles.shape[0]
    mass = f.root.mass[()]
    p_ref = f.root.pz[()]
    sn = f.root.s_n[()] #period length
    tn = f.root.tlen[()] #cumulative tracked length for this file

    f.close()
    
    #ID = np.arange(5000)
    
    header = dict()
    header['n_part'] = npart
    header['mass'] = mass
    header['p_ref'] = p_ref
    header['s_val'] = sn
    header['t_len'] = tn
    
    part_vals = []
    
    #separate lost particles
    for particle in particles:
        val = particle[6]
        if val in ID:
            part_vals.append(particle)
            
    return np.asarray(part_vals)

def get_particles(inputfile, lost=None, lostlist=None):
    
    '''
    Reads an input file and returns a numpy array of particles and a dictionary of root values. 
    If lost particles specified, then returns those separately.
    
    Arguments:
        inputfile (str): path to a .h5 file containing particle diagnostics.
        lost (Optional[bool]): True if some particles have been lost. Defaults to None.
        lostlist (Optional[list]): List of particle IDs that were lost during the simulation. Defaults to None.
        
    Returns:
        header (dict): dictionary of meta data regarding the particle array
        particles (ndArray): array of particle data [x, x', y, y', cdt, dp, ID] for each particle    
    
    '''
    
    f = tables.openFile(inputfile, 'r')
    particles = f.root.particles.read()
    
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
    
    if lost:
    
        #define new lists
        lost_particles = []
        kept_particles = []
        lost_counter = lostlist
    
        #separate lost particles
        for index,particle in enumerate(particles):
            val = particle[6]
            if val in lostlist:
                lost_particles.append(particle)
            else:
                kept_particles.append(particle)
    
        return header, np.asarray(kept_particles), np.asarray(lost_particles)
    
    else:
        
        return header, particles
    

def get_file_list(opts, prefix='particles'):

    '''
    Returns a list of files of the form 'particles*.h5' from the current directory 

    Arguments:
        opts (options.Options): Synergia options instance specifying relative path
        prefix (Optional[str]): Used to distinguish different .h5 file prefixes. Default is 'particles'
    
    '''
    
    if not opts.relpath:
        #If no relative path specified, check current directory 
        files = os.listdir(os.getcwd())
    else:
        #Otherwise check specified relative path
        path = os.path.join(os.getcwd(),opts.relpath)
        files = [os.path.join(path,fn) for fn in os.listdir(path)]

    pfiles = []

    #basic filtering for particles*.h5 files
    for filename in files:
            if filename.find(str(prefix)) > -1 and filename.endswith('.h5'):
                pfiles.append(filename)
    
    return np.sort(pfiles)
    
def get_lost_particle_list(opts):
    '''
    Returns a list of particle IDs corresponding to lost particles from inspection of the output files
    
    Arguments:
        opts (options.Options): Synergia options instance specifying relative path
        
    Returns:
        lost_vals (ndArray): array of particle IDs corresponding to particles that were lost
    
    '''
    files = get_file_list(opts)
    
    #compare first file output to final file output
    header1, particles1 = get_particles(files[0])
    
    header2, particles2 = get_particles(files[-1])
    
    if header1['n_part'] == header2['n_part']:
        #we have no lost particles
        lost_vals = []
    
    else:
        indices1 = particles1[:,6]
        indices2 = particles2[:,6]
    
        combined_index = np.append(indices1,indices2)
        s = np.sort(combined_index)
        ci_shared = s[s[1:] == s[:-1]] #list of those that remain
        ci_full = [int(ind) for ind in np.unique(combined_index)] #full list
        lost_vals = np.delete(ci_full, ci_shared) #lost values

        if not len(lost_vals) == len(ci_full) - len(ci_shared):
            #Print a sequence of warnings and some debug information
            print "Warning: Length of lost list is not equal to number of lost particles!"
            print "{} values are shared out of {} total values.".format(len(ci_shared),len(ci_full))
            print "Therefore there are {} lost values.".format(len(ci_full)-len(ci_shared))
            print "However I caclulate the length of the lost array to be {}.".format(len(lost_vals))
    
    return lost_vals
    
    
def get_twiss(lattice_simulator):
    '''
    Returns an array of twiss parameters versus longitudinal coordinate 's' for a given lattice.
    
    Arguments:
        lattice_simulator (synergia.simulation.lattice_simulator) - a Synergia lattice simulator
    
    Returns:
        twiss (ndarray): Array of twiss values with structure [s,betax,alphax,gammax,betay,alphay,gammay]
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
    
    This represents an improvement in resolution of the lattice, as latticefunctions are calculated slice by slice. 
    See lfplot.get_sliced_lattice_functions(), as it represents the same operation. While this in many ways deprecates
    the original "get_twiss()" parameterization, its relevant for specific instances where you want a "slice-independent"
    view of the lattice functions to compare against other possible implementations.
    
    Arguments:
        lattice_simulator (synergia.simulation.lattice_simulator) - a Synergia lattice simulator
    
    Returns:
        twiss (ndarray): Array of twiss values with structure [s,betax,alphax,gammax,betay,alphay,gammay]
    '''
    
    lattice = lattice_simulator.get_lattice()
    
    twiss = []
    
    for aslice in lattice_simulator.get_slices():
        temp = []
        
        lattice_functions=lattice_simulator.get_lattice_functions(aslice)
        
        #if no twiss thus far, just add the first entry
        if twiss == []:
            gamma_x = (1 + lattice_functions.alpha_x**2)/lattice_functions.beta_x
            gamma_y = (1 + lattice_functions.alpha_y**2)/lattice_functions.beta_y
            temp = [lattice_functions.arc_length, lattice_functions.beta_x, lattice_functions.alpha_x, gamma_x, lattice_functions.beta_y, lattice_functions.alpha_y, gamma_y]
            twiss.append(temp)            
    
        #if many slices in, check that the next slice is at a new s position before appending twiss parameters
        elif  not twiss[-1][0] == lattice_functions.arc_length:
            
            gamma_x = (1 + lattice_functions.alpha_x**2)/lattice_functions.beta_x
            gamma_y = (1 + lattice_functions.alpha_y**2)/lattice_functions.beta_y
            temp = [lattice_functions.arc_length, lattice_functions.beta_x, lattice_functions.alpha_x, gamma_x, lattice_functions.beta_y, lattice_functions.alpha_y, gamma_y]
            twiss.append(temp)
        else:
            pass
    
    #at the end, pre-pend final twiss values @ beginning with s-value 0
    #twiss_init = twiss[-1]
    #twiss_init[0]=0.0
    #twiss.insert(0,twiss_init)
    
    #verify that we haven't inserted an improper 0 value at the end
    if twiss[-1][0] < twiss[-2][0]:
            twiss = twiss[:-1][:]
        
    return np.asarray(twiss)
    

def get_human_coords(filelist, lost=None,lostlist=None, plotlost=False, num=None):
    '''
    Returns a numpy array of non-normalized (human) coordinate vectors from .h5 files in filelist.
    
    Arguments:
        filelist (str):  A list of .h5 files containing particle array information
        lost (Optional[bool]): True if some particles have been lost. Defaults to None.
        lostlist (Optional[list]): List of particle IDs that were lost during the simulation. Defaults to None.
        plotlost (Optional[bool]): True if the lost particle coordinates are desired. Default to False.
        num (Optional[int]): If specified, only returns the first num particles. Default to None (returns all).
    
    Outputs:
        full_coords (ndarray): A numpy array with dimension #turns x #particles x #transverse coordinates(4).
    
    
    Note that the command norm[B][A] returns vector of coordinates for particle A at turn B.
    
    '''
    full_coords = [] #full_coords is a list of arrays of particle coordinates
    if lost:
        max_num = len(lostlist) #maximum # of particles in case of particle lost
        #print max_num
    
    for index,fileName in enumerate(filelist):
        inputfile = fileName         
        
        if lost:
            header, particles, lost_particles = get_particles(inputfile, lost,lostlist)
        
            x = lost_particles[:,coords['x']] #units m
            xp = lost_particles[:,coords['xp']] #unitless
            y = lost_particles[:,coords['y']] #units m
            yp = lost_particles[:,coords['yp']] #unitless 
            #stack arrays then transpose to get array of coordinate vectors
            pID =  lost_particles[:,6]
            lost_particles_cut = np.vstack((x,xp,y,yp)).T
        
        else:
            header, particles = get_particles(inputfile)
        
        x = particles[:,coords['x']] #units m
        xp = particles[:,coords['xp']] #unitless
        y = particles[:,coords['y']] #units m
        yp = particles[:,coords['yp']] #unitless
        pID =  particles[:,6]
        #stack arrays then transpose to get array of coordinate vectors
        particles_cut = np.vstack((x,xp,y,yp)).T
        
        
        if plotlost:
            if lost_particles.size == 0:
                coord = np.zeros((max_num,4)) #append all zeros
                pass #don't update normalized coordinates because lost particle array has no more particles!
            else:
                partial_coords = lost_particles_cut #get partial normal coordinates
                 
                num_left = partial_coords.shape[0]
                #num_left = header['n_part']
                num_0 = max_num - num_left
                if num_0 > 0:
                    nzeros = np.zeros((num_0,4)) #construct an array of zeros to append to partial_coords
                    coord = np.vstack((partial_coords,nzeros))
                else:
                    coord = partial_coords #we have all particles
        else:
            coord = particles_cut
        
        #only append first num of tracked particles if not plotting all
        if num:
            full_coords.append(coord[:num])
        else:  
            full_coords.append(coord) #added flatten here

    return np.asarray(full_coords)
    
def get_normalized_coords(filelist, twiss, lost=None,lostlist=None, plotlost=False, num=None, ID=None):
    
    '''
    Returns a numpy array of normalized coordinate vectors obtained from .h5 files in filelist
    
    Arguments:
        filelist (str):  A list of .h5 files containing particle array information
        twiss (array-like): An array of twiss parameters
        lost (Optional[bool]): True if some particles have been lost. Defaults to None.
        lostlist (Optional[list]): List of particle IDs that were lost during the simulation. Defaults to None.
        plotlost (Optional[bool]): True if the lost particle coordinates are desired. Defaults to False.
        num (Optional[int]): If specified, only returns the first num particles. Defaults to None (returns all).
        ID (Optional[int]): If specified, only returns coordinates for particle with ID. Defaults to None.
    
    Outputs:
        norms (ndarray): A numpy array with dimension #turns x #particles x #transverse coordinates(4).
    
    Returns a numpy array with dimensions #turns x #particles x #transverse coordinates(4).
    norms[B][A] returns vector of coordinates for particle A at turn B.
    '''
    
    norms = [] #norms is a list of arrays of particle coordinates
    if lost:
        max_num = len(lostlist) #maximum # of particles in case of particle lost
        #print max_num
    
    for index,fileName in enumerate(filelist):
        inputfile = fileName         
        
        if lost:
            header, particles, lost_particles = get_particles(inputfile, lost,lostlist)
        else:
            header, particles = get_particles(inputfile)
        
        if plotlost:
            if lost_particles.size == 0:
                norm_coords = np.zeros((max_num,4)) #append all zeros
                pass #don't update normalized coordinates because lost particle array has no more particles!
            else:
                partial_coords = normalized_coordinates(header,lost_particles,twiss) #get partial normal coordinates 
                num_left = partial_coords.shape[0]
                #num_left = header['n_part']
                num_0 = max_num - num_left
                if num_0 > 0:
                    nzeros = np.zeros((num_0,4)) #construct an array of zeros to append to partial_coords
                    norm_coords = np.vstack((partial_coords,nzeros))
                else:
                    norm_coords = partial_coords #we have all particles
        else:
            norm_coords = normalized_coordinates(header, particles, twiss)
        
        #only append first num of tracked particles if not plotting all
        if num:
            norms.append(norm_coords[:num])
        else:  
            norms.append(norm_coords) #added flatten here

    return np.asarray(norms)



def get_single_particle_elliptic_invariants(filelist, twiss, opts, lost, num=1):
    '''
    Returns an array of single particle invariants (Hamiltonian with elliptic potential)
    
    Arguments:
        filelist (str):  A list of .h5 files containing particle array information
        twiss (array-like): An array of twiss parameters
        lost (Optional[bool]): True if some particles have been lost. Defaults to None.
        lostlist (Optional[list]): List of particle IDs that were lost during the simulation. Defaults to None.
        num (Optional[int]): If num=2, return 2nd invariant. Otherwise return 1st invariant. Defaults to 1.
        
    Returns:
        invariant (ndarray): Array of specified invariant value calculated for each file (e.g. each turn)
    
    '''

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

    return np.asarray(invariant)

def get_adjusted_elliptic_invariants(filelist, twiss, opts, lost, num=1):

    '''
    
    Returns an array of single particle invariants (Hamiltonian with elliptic potential)
    adjusted for off momentum particles with a properly chromaticity corrected lattice.
    
    See Webb et al. "Chromatic and Dispersive Effects in Nonlinear Integrable Optics"
    for a discussion of the requirements for the transformation of nonlinear strength parameter
    t to be applied in the calculation of the adjusted invariant.
    
    Arguments:
        filelist (str):  A list of .h5 files containing particle array information
        twiss (array-like): An array of twiss parameters
        lost (Optional[bool]): True if some particles have been lost. Defaults to None.
        lostlist (Optional[list]): List of particle IDs that were lost during the simulation. Defaults to None.
        num (Optional[int]): If num=2, return 2nd invariant. Otherwise return 1st invariant. Defaults to 1.
        
    Returns:
        invariant (ndarray): Array of specified invariant value calculated for each file (e.g. each turn)
    
    '''

    invariant = [] #invariant is a list of arrays of macroparticles
    
    for index,fileName in enumerate(filelist):
        inputfile = fileName
        
        if lost:
            header, particles, lost_particles = get_particles(inputfile, lost)
        else:
            header, particles = get_particles(inputfile, lost)
        
        #normalize coordinates    
        norm_coords = normalized_coordinates(header, particles, twiss)
        
        #grab deltap/p array
        delta_vals = particles[:,5]
        
        t_array = modified_T_array(opts.t,delta_vals,opts)
        base_t = opts.t*np.ones(len(delta_vals))
        t_ratio = np.divide(t_array,base_t)
        
        #construct elliptic coordinates for that file
        u,v = elliptic_coordinates(norm_coords, opts)
        
        
        if num == 1:
            #calculate the first invariant
            vals = single_particle_hamiltonian(norm_coords) + np.multiply(t_ratio,elliptic_hamiltonian(u,v,opts))         
        elif num == 2:
            #calculate the second invariant
            vals = second_invariant(norm_coords, u,v, opts)
        else:
            #otherwise calculate the first invariant
            print "Improper invariant number specified. Calculating 1st invariant by default."
            vals = single_particle_hamiltonian(norm_coords) + np.multiply(t_ratio,elliptic_hamiltonian(u,v,opts))
        
        invariant.append(vals)
        
    return np.asarray(invariant)

def get_toy_twiss(opts):
    '''
    Returns a fixed twiss array for all lattice 's' values.
    
    For an idealized nonlinear magnet, the beta and alpha values can be specified given a phase advance
    and length of the underlying drift. This script assumes that such values were comptued and given to
    an Options instance for computation of normalized coordinates with an idealized lattice.
    
    Arguments:
        opts (options.Options): Synergia options instance specifying the necessary twiss parameters.
    
    Returns:
        twiss (ndarray): Array of twiss values with form [s, betax, alphax, gammax, betay, alphay, gammay]
    
    '''

    
    svals = np.linspace(0,opts.lattice.get_length(),100)
    
    length = len(svals)
    twiss = np.zeros((length,7))    
    
    
    #opts.betae = 1.*opts.betae
    gamma = (1 + opts.alphae**2)/(opts.betae)
    
    
    twiss[::,0] = svals
    twiss[::,1] = [opts.betae for i in range(length)]
    twiss[::,2] = [opts.alphae for i in range(length)]
    twiss[::,3] = [gamma for i in range(length)] 
    twiss[::,4] = [opts.betae for i in range(length)]
    twiss[::,5] = [opts.alphae for i in range(length)]
    twiss[::,6] = [gamma for i in range(length)] 

    return twiss



################################## COORDINATE TRANSFORMS ####################################
def modified_T(t, delta, opts):
    """
    Computes the modified nonlinear strength for use with Chromaticity correction.
    
    See Webb et al. "Chromatic and Dispersive Effects in Nonlinear Integrable Optics"
    for a discussion of the requirements for the transformation of nonlinear strength parameter
    t to be applied in the calculation of the adjusted invariant.
    
    Arguments:
        t (float): The strength parameter for the nonlinear magnet
        delta (float): The relative momentum deviation of the particle.
        opts (options.Options): A Synergia options instance specifying the needed phase advance quanties
    
    Returns:
        t/correction (float): A corrected effective strength for the nonlinear magnet
    
    """
    
    mu0 = opts.full_tune
    nu0 = opts.tune
    Ch = opts.Ch
    correction = 1. - (mu0*delta*Ch/nu0)
    
    return t/correction
def modified_T_array(t, delta_array, opts):
    
    """
    Computes the modified nonlinear strength for use with Chromaticity correction for an array of particles.
    
    See Webb et al. "Chromatic and Dispersive Effects in Nonlinear Integrable Optics"
    for a discussion of the requirements for the transformation of nonlinear strength parameter
    t to be applied in the calculation of the adjusted invariant.
    
    Arguments:
        t (float): The strength parameter for the nonlinear magnet
        delta_array (ndarray): An array containing the relative momentum deviations for a group of particles.
        opts (options.Options): A Synergia options instance specifying the needed phase advance quanties
    
    Returns:
        t/correction (ndarray): Array of corrected effective strength values for each particle
    
    """
    
    mu0 = opts.full_tune
    nu0 = opts.tune
    Ch = opts.Ch
    
    t_array = t*np.ones(len(delta_array))
    
    correction = 1.*np.ones(len(delta_array)) - (mu0*delta_array*Ch/nu0)
    
    return np.divide(t_array,correction)
    
def phase_advance(p1, p2):
    '''
    Returns the angle between two vectors.
    
    This script is used to compute an effective phase advance by comparing an initial and final position in phase space.
    
    Arguments:
        p1 (ndarray): Array containing phase space coordinates (e.g. [x, x'])
        p2 (ndarray): Array containing phase space coordinates (e.g. [x, x'])
        
    Returns:
        guess_angle (float): The angle between the input vectors, given in radians. Assumes clockwise rotation.
    
    '''
    
    norm1 = np.linalg.norm(p1)
    norm2 = np.linalg.norm(p2)
    
    #Determine the quadrant - the integer corresponds to the quadrant
    q1 = 1*(p1[0] > 0) & (p1[1] > 0) + 2*(p1[0] < 0) & (p1[1] > 0) + 3*(p1[0] < 0) & (p1[1] < 0) + 4*(p1[0] > 0) & (p1[1] < 0)
    q2 = 1*(p2[0] > 0) & (p2[1] > 0) + 2*(p2[0] < 0) & (p2[1] > 0) + 3*(p2[0] < 0) & (p2[1] < 0) + 4*(p2[0] > 0) & (p2[1] < 0)

    #calculate dot product
    product = np.dot(p1,p2)
    guess_angle = np.arccos(product/(norm1*norm2))
    
    if not q1==q2:
        #quadrants change
        if q2 > q1 or (q2==1 and q1 ==4):
            #we have clockwise rotation across quadrants
            return guess_angle
        else:
            #we have counter clockwise rotation across quadrants
            return guess_angle
            #return (2*np.pi* - guess_angle)
    else:
        #same quadrant, so check x-value
        #clockwise guess
        guess_pos_p2 = (norm2/norm1)*(p1[0]*np.cos(guess_angle) - p1[1]*np.sin(guess_angle)) 
        #counterclockwise guess
        guess_neg_p2 = (norm2/norm1)*(p1[0]*np.cos(guess_angle) + p1[1]*np.sin(guess_angle))
        #print guess
        if np.abs(guess_pos_p2 - p2[0]) < np.abs(guess_neg_p2 - p2[0]):
            #positive guess is closer, clockwise
            return guess_angle
        else:
            #negative closer, so counterclockwise
            return guess_angle
    
    
def normalized_coordinates(header, particles, twiss, offset=None, ID=None):
    '''
    Return an array of particles (fixed s) in linear normalized transverse coordinates.
    
    This function computes the normalized coordinates given an array of particle trace space
    coordinates and twiss parameters. Particle coordinates are transformed in the following way
    
        Input Coordinates - x, x', y, y'
        Output Coordinates - sqrt(betax)*x, 
                            (1/sqrt(betax))*(betax*x' + alpha*x), 
                            sqrt(betay)*y, 
                            (1/sqrt(betay))*(betay*y' + alpha*y)
                        
    Arguments:
        header (str):  A list of .h5 files containing particle array information
        particles (ndarray): An array of particle trace space coordinates
        twiss (array-like): An array of twiss parameters
        offset (Optional[float]): Specify an svalue offset to be applied to all coordinates. Defaults to None.
        ID (Optional[int]): If specified, only returns coordinates for particle with ID. Defaults to None.
    
    Outputs:
        norms (ndarray): A numpy array with dimension #turns x #particles x #transverse coordinates(4).
        
    
    '''
    
    sval = header['s_val']
    
    #Must be careful to use the appropriate s_val, or else normalizaztion will be off
    #Thus we force periodicity!
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
     
    #Keep in mind that alpha is flipped at the end of the NL element
    alphax = 1*alphax
    alphay = 1*alphay
    
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
    


def elliptic_coordinates(normalized, opts):
    '''
    Return the an array of elliptic coordinate values u and v for a list of particles.
    
    This function computes the elliptic coordinates given an array of particle trace space
    coordinates and twiss parameters. Particle coordinates are transformed in the following way
        Input Coordinates - normal coordinates x, px, y, py
        Output Coordinates - u, v - vectors of length equal to the # of particles
    
    Arguments:
        normalized (array-like): Array of particles normalized coordinates
        opts (options.Options): A Synergia options instance specifying nonlinear magnet parameters
        
    Returns:
        [u,v] (ndarray): tuple containing arrays of elliptic coordinates  
    '''

    c = opts.c 
        
    x = normalized[:,0]
    y = normalized[:,2]
    
    #first need to adjust x and y by the c factor
    x = x*1.0/c
    y = y*1.0/c
    
    u = 0.5*(np.sqrt((x + 1.)**2 + y**2) + np.sqrt((x -1.)**2 + y**2))
    v = 0.5*(np.sqrt((x + 1.)**2 + y**2) - np.sqrt((x -1.)**2 + y**2))
    
    return [u,v]
    
    
def elliptic_hamiltonian(u,v, opts):
    '''
    Returns arrays of values for the first elliptic invariant (Hamiltonian) for a system with the IOTA nonlinear potential.
    
    Arguments:
        u (ndarray): Array of particles first elliptic coordinate
        v (ndarray): Array of particles second elliptic coordinate
        opts (options.Options): A Synergia options instance specifying nonlinear magnet parameters
        
    Returns:
        kfac*elliptic (ndarray): Array of values for the elliptic Hamiltonian for each particle  
    '''

    
    t = opts.t 
    c = opts.c
    
    f2u = u * np.sqrt(u**2 -1.) * np.arccosh(u)
    g2v = v * np.sqrt(1. - v**2) * (-np.pi/2 + np.arccos(v))
    
    elliptic = (f2u + g2v) / (u**2 - v**2)
    kfac = -1.*t*c*c
    
    return kfac*elliptic
    
    
def second_invariant(normalized, u,v, opts):
    '''
    Returns arrays of values for the second elliptic invariant (Hamiltonian) for a system with the IOTA nonlinear potential.
    
    Arguments:
        normalized (array-like): Array of particles normalized coordinates (x, px, y, py)
        u (ndarray): Array of particles first elliptic coordinate
        v (ndarray): Array of particles second elliptic coordinate
        opts (options.Options): A Synergia options instance specifying nonlinear magnet parameters
        
    Returns:
        invariant (ndarray): Array of values for the second invariant for each particle  
    '''
    
    t = -1.*opts.t
    c = 1.*opts.c
    
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
    
    invariant = (p_ang + p_lin) + 2.*(c**2) * (fu * v**2 + gv * u**2)/(u**2 - v**2)
    
    return invariant
       

def single_particle_invariant(header, particles, twiss):
    '''
    Computes the Courant-Synder invariant for an array of particles.
    
    Arguments:
        header (str):  A list of .h5 files containing particle array information
        particles (ndarray): An array of particle trace space coordinates
        twiss (array-like): An array of twiss parameters
        
    Returns:
        inv2 (ndarray): Array of values for the Courant-Synder invariant
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
    
    inv2 = invariantx + invarianty
    
    return inv2
    
def get_invariants(filelist, twiss, lost):
    '''
    Computes the Courant-Synder invariant for a series of arrays of particles from a list of files.
    
    Arguments:
        filelist (str):  A list of .h5 files containing particle array information
        twiss (array-like): An array of twiss parameters
        lost (Optional[bool]): True if some particles have been lost. Defaults to None.
        
    Returns:
        invariant (ndarray): Array of invariant values
    '''

    invariant = [] #hamiltonian is a list of arrays of macroparticles
    
    for index,fileName in enumerate(filelist):
        inputfile = fileName
        
        if lost:
            header, particles, lost_particles = get_particles(inputfile, lost)
        else:
            header, particles = get_particles(inputfile, lost)
            
        norm_coords = normalized_coordinates(header, particles, twiss)
        invariant.append(single_particle_hamiltonian(norm_coords))
        
    return np.asarray(invariant)


def single_particle_hamiltonian(normalized, ID=None):
    '''
    Computes the single particle Hamiltonian for an array of particles in absence of nonlinearities.
    
    Arguments:
        normalized (array-like): Array of particles normalized coordinates (x, px, y, py)
        ID(Optional[int]): Specifies a single particle for computation. Defaults to None.
        
    Returns:
        quadradtic (array-like): Array of Hamiltonian values for specified particles
    
    '''
    
    x = normalized[:,0]
    px = normalized[:,1]
    y = normalized[:,2]
    py = normalized[:,3]
    
    
    if ID:
        #quadratic is the basis for calculating the hamiltonian, in absence of nonlinear couplings
        quadratic = 0.5* (px[ID]**2 + py[ID]**2) + 0.5*(x[ID]**2 + y[ID]**2)
    else:
        quadratic = 0.5* (px**2 + py**2) + 0.5*(x**2 + y**2)
        
    return quadratic
    

def get_hamiltonians(filelist):
    
    '''
    
    Computes the Hamiltonian from particle coordinates obtained from .h5 files
    
    Arguments:
        filelist (list): A list of .h5 files containing particle array information
    
    '''

    hamiltonian = [] #hamiltonian is a list of arrays of macroparticles
    
    for index,fileName in enumerate(filelist):
        inputfile = fileName
        header, particles = get_particles(inputfile)
        hamiltonian.append(single_particle_hamiltonian(header,particles))
        
    return np.asarray(hamiltonian)


#################################Called Scripts#################################################################

def single_turn_phase_advance(files, twiss, dim='x', nParticles=1000, indices=[0,1]):
    '''
    Compute the phase advance, modulo 2 pi, for particles between their outputs to two separate files.
    
    This function compiles normalized coordinates and returns the phase advance between those
    particles from the first file specified by the turns parameter to the second file specied.
    
    Arguments:
        filelist (str):  A list of .h5 files containing particle array information
        twiss (array-like): An array of twiss parameters
        dim (Optional[str]): Specifies plane ('x' or 'y') for computing phase advance. Defaults to 'x'.
        nParticles (Optional[int]): Specifies the number of particles to include in the return array. Defaults to [0,999).
        indices (Optional[int]): Specifies the start and end file to look at. Defaults to [0,1]
    
    
    Returns:
        Phases (array): Array of phase advances, nParticles in length
    
    '''
    
    phases = []
    
    needed_files = [files[ind] for ind in indices] #slice complete file list by indices
    
    norm_coords = get_normalized_coords(files,twiss)
    
    for ind in range(nParticles):
        if dim == 'x':
            p1 = norm_coords[indices[0]][ind,(0,1)]
            p2 = norm_coords[indices[1]][ind,(0,1)]
            phases.append(phase_advance(p1,p2)/(2.*np.pi))
        if dim == 'y':
            p1 = norm_coords[indices[0]][ind,(2,3)]
            p2 = norm_coords[indices[1]][ind,(2,3)]
            phases.append(phase_advance(p1,p2)/(2.*np.pi))  
    
    
    return phases
    
def plot_Poincare(opts, noTwiss=False):
    '''
    Plot a Poincare section in the desired normalized coordinates using the full Twiss of the lattice.
    
    This function extrapolates twiss parameters from a lattice simulator and grabs diagnostic files from
    an output directory in order to compute and plot normalized coordinates (e.g. x-px or xy plots).
    
    Arguments:
        opts (options.Options):  A Synergia options instance containing plotting parameters and path information
    
    '''
    
    opts.hcoord = opts.plots[0]
    opts.vcoord = opts.plots[1]
    
    files = get_file_list(opts)
    lost = get_lost_particle_list(opts)

    if len(lost) > 0:
        #we have lost particles
        opts.lost = lost #store these in opts.lost
        lost = True #make lost a simple flag
    
    twiss = get_sliced_twiss(opts.lattice_simulator)
    
    if opts.plot_lost:
         pArray = get_normalized_coords(files,twiss,lost,opts.lost,True)
        
    else:
        pArray = get_normalized_coords(files,twiss,lost,opts.lost)
    
    plot_P(pArray, opts) 
    
def toy_plot_Poincare(opts):
    '''
    Plot a Poincare section in the desired normalized coordinates using fixed externally provided Twiss parameters.
    
    This function uses pre-determined "toy" Twiss parameters from get_toy_twiss() and grabs diagnostic files from
    an output directory in order to compute and plot normalized coordinates (e.g. x-px or xy plots).
    
    Arguments:
        opts (options.Options):  A Synergia options instance containing plotting parameters and path information.
                                 Can also specify opts.plot_lost to look at trajectories that are lost during the
                                 course of the simulation.
    
    '''
    
    opts.hcoord = opts.plots[0]
    opts.vcoord = opts.plots[1]
    
    files = get_file_list(opts)
    lost = get_lost_particle_list(opts)
    
    if len(lost) > 0:
        #we have lost particles
        opts.lost = lost #store these in opts.lost
        lost = True #make lost a simple flag

    twiss = get_toy_twiss(opts)
    
    
    if opts.plot_lost:
         pArray = get_normalized_coords(files,twiss,lost,opts.lost,True)
        
    else:
        pArray = get_normalized_coords(files,twiss,lost,opts.lost)
    
    plot_P(pArray, opts) 


def plot_invariant(opts):
    '''
    Plot the Courant-Synder invariant over a specified number of turns and number of particles
    
    This function extrapolates twiss parameters from a lattice simulator and grabs diagnostic files from
    an output directory in order to compute normalized coordinates and then the Courant-Synder Invariant.
    
    Arguments:
        opts (options.Options):  A Synergia options instance containing plotting parameters and path information.
                                 This object also points to the lattice_simulator and other simulation specific details.
    
    '''
    
    opts.hcoord = 'turn #'
    opts.vcoord = 'J(p,q)'
    opts.elliptic = False
    
    files = get_file_list(opts)
    twiss = get_sliced_twiss(opts.lattice_simulator)
    lost = get_lost_particle_list(opts)
    jArray = get_invariants(files, twiss, lost)
    plot_J(jArray,opts)
    
def plot_elliptic_invariant(opts):
    '''
    Plot an elliptic invariant over a specified number of turns and number of particles.
    
    This function extrapolates twiss parameters from a lattice simulator and grabs diagnostic files from
    an output directory in order to compute normalized coordinates and then one of the two elliptic invariants.
    The opts.num parameter is used to specify which invariant is plotted.
    
    Arguments:
        opts (options.Options):  A Synergia options instance containing plotting parameters and path information.
                                 This object also points to the lattice_simulator and other simulation specific details.
                                 In particular, use opts.num=2 to plot the second invariant.
    '''
    
    opts.hcoord = 'turn #'
    opts.elliptic = True
    
    if opts.num:
        num = opts.num
    else:
        num = 1 #set as default
        opts.num = 1
        
    if num == 2:
        opts.vcoord = 'I(p,q)'
    else:
        opts.vcoord = 'H(p,q)' #First invariant (H) is the default!
    
    files = get_file_list(opts)
    twiss = get_sliced_twiss(opts.lattice_simulator)
    lost = get_lost_particle_list(opts)
    jArray = get_single_particle_elliptic_invariants(files, twiss, opts, lost, num)
    plot_J(jArray,opts)
    
    

def toy_plot_elliptic_invariant(opts):
    '''
    Plot an elliptic invariant over a specified number of turns and number of particles.
    
    Unlike "plot_elliptic_Invariant", this file grabs a fixed twiss value from get_toy_twiss! It then grabs 
    diagnostic files from an output directory in order to compute normalized coordinates and then one of the 
    two elliptic invariants. The opts.num parameter is used to specify which invariant is plotted.
    
    Arguments:
        opts (options.Options):  A Synergia options instance containing plotting parameters and path information.
                                 This object also points to the lattice_simulator and other simulation specific details.
                                 In particular, use opts.num=2 to plot the second invariant.
    
    '''
    
    opts.hcoord = 'turn #'
    opts.elliptic = True
    
    if opts.num:
        num = opts.num
    else:
        num = 1 #set as default
        opts.num = 1
        
    if num == 2:
        opts.vcoord = 'I(p,q)'
    else:
        opts.vcoord = 'H(p,q)' #First invariant (H) is the default!
    
    files = get_file_list(opts)
    twiss = get_toy_twiss(opts)
    lost = get_lost_particle_list(opts)
    jArray = get_single_particle_elliptic_invariants(files, twiss, opts, lost, num)
    plot_J(jArray,opts)
    

def toy_calc_bunch_H(bunch, opts, elliptic = True):
    '''
    Calculate the invariants for a group of particles.
        
    This script accepts as input both ndarrays of bunches, and Synergia bunches.
    Unlike "calc_bunch_H", this file grabs a fixed twiss value from get_toy_twiss in 
    order to compute normalized coordinates and then the invariants.
    
    Arguments:
        bunch (synergia.bunch.bunch.Bunch/ndarray) : A Synergia bunch object or an array of particles - [x, x', y, y', cdt, dp, ID]
        opts (options.Options):  A Synergia options instance containing plotting parameters and path information.
                                 This object also points to the lattice_simulator and other simulation specific details.   
        elliptic (Optional[bool]): If True, compute elliptic invariants. Defaults to True.
    
    Returns:
        hArray, iArray (ndarray): Arrays of invariant values - if elliptic=False, iArray defaults to zeroes.
    '''
    #We'd like to use this for test bunches as well as synergia bunches
    if type(bunch) == synergia.bunch.bunch.Bunch:
        particles = bunch.get_local_particles()
    elif type(bunch) == np.ndarray:
        particles = bunch
    
    twiss = get_toy_twiss(opts)
    header= {}
    header['s_val'] = 0.
    norm_coords = normalized_coordinates(header, particles, twiss)
    
    if elliptic:
        #calculate the elliptic quantities
        u,v = elliptic_coordinates(norm_coords, opts)
        hArray = single_particle_hamiltonian(norm_coords) + elliptic_hamiltonian(u,v,opts)  
        iArray = second_invariant(norm_coords,u,v,opts)
        return hArray, iArray
        
    else:
        #calculate the regular Hamiltonian
        hArray = single_particle_hamiltonian(norm_coords)
        iArray = np.zeros(len(particles))
        return hArray, iArray
        

def toy_calc_H_and_ID(bunch, opts, elliptic = True):
    '''
    Calculate the first invariant for a group of particles and returns a corresponding particle ID list.
        
    This script accepts as input both ndarrays of bunches, and Synergia bunches.
    Unlike "calc_bunch_H", this file grabs a fixed twiss value from get_toy_twiss in 
    order to compute normalized coordinates and then the invariants.
    
    Arguments:
        bunch (synergia.bunch.bunch.Bunch/ndarray) : A Synergia bunch object or an array of particles - [x, x', y, y', cdt, dp, ID]
        opts (options.Options):  A Synergia options instance containing plotting parameters and path information.
                                 This object also points to the lattice_simulator and other simulation specific details.
        elliptic (Optional[bool]): If True, compute elliptic invariants. Defaults to True.
    
    Returns:
        hID (ndarray): arrays of first invariant values and particle ID values, respectively
    '''
    
    #We'd like to use this for test bunches as well as synergia bunches
    if type(bunch) == synergia.bunch.bunch.Bunch:
        particles = bunch.get_local_particles()
    elif type(bunch) == np.ndarray:
        particles = bunch
    
    twiss = get_toy_twiss(opts)
    header= {}
    header['s_val'] = 0.
    norm_coords = normalized_coordinates(header, particles, twiss)
    
    ID_vals = particles[:,6]
    IDV = ID_vals.reshape(len(ID_vals),1) #force a reshape for stacking purposes
    
    
    if elliptic:
        u,v = elliptic_coordinates(norm_coords, opts)
        hArray = single_particle_hamiltonian(norm_coords) + elliptic_hamiltonian(u,v,opts)  
        iArray = second_invariant(norm_coords,u,v,opts)
        
        hID = np.transpose(np.vstack((ID_vals,hArray)))
        iID = np.transpose(np.vstack((ID_vals,iArray)))
        
        return hID, iID
        
    else:
        #calculate the regular Hamiltonian
        hArray = single_particle_hamiltonian(norm_coords)
        iArray = np.zeros(len(particles))
        
        hID = np.transpose(np.vstack((ID_vals,hArray)))
        iID = np.transpose(np.vstack((ID_vals,iArray)))        

        return hID, iID

def calc_bunch_H(bunch, opts,header = None, elliptic = True):
    '''
    Calculate the invariants for a group of particles for any header.
        
    This script accepts as input both ndarrays of bunches, and Synergia bunches. It grabs the full twiss of the lattice
    from the provided opts.lattice_simulator in order to compute normalized coordinates and then the invariants.
    
    Arguments:
        bunch (synergia.bunch.bunch.Bunch/ndarray) : A Synergia bunch object or an array of particles - [x, x', y, y', cdt, dp, ID]
        opts (options.Options):  A Synergia options instance containing plotting parameters and path information.
                                 This object also points to the lattice_simulator and other simulation specific details.
        header (dict): A Python dictionary containing meta-data for a particle diagnostics file. Specifies 's_val' value.
                        Defaults to a dummy header with 's_val' = 0.
        elliptic (Optional[bool]): If True, compute elliptic invariants. Defaults to True.
    
    Returns:
        hArray, iArray (ndarray): Arrays of invariant values - if elliptic=False, iArray defaults to zeroes.
    '''
    
    #We'd like to use this for test bunches as well as synergia bunches
    if type(bunch) == synergia.bunch.bunch.Bunch:
        particles = bunch.get_local_particles()
    elif type(bunch) == np.ndarray:
        particles = bunch
    
    twiss = get_sliced_twiss(opts.lattice_simulator)
    
    #if no header provided, create dummy header with s_val = 0!
    if not header:
        header = {}
        header['s_val'] = 0
    
    norm_coords = normalized_coordinates(header, particles, twiss)
    
    if elliptic:
        u,v = elliptic_coordinates(norm_coords, opts)
        hArray = single_particle_hamiltonian(norm_coords) + elliptic_hamiltonian(u,v,opts)  
        iArray = second_invariant(norm_coords,u,v,opts)
        return hArray, iArray
        
    else:
        #calculate the regular Hamiltonian
        hArray = single_particle_hamiltonian(norm_coords)
        iArray = np.zeros(len(particles))
        return hArray, iArray

def calc_elliptic_Invariant(opts, elliptic=True):
    '''
    Returns both invariants for NL elliptic potential over a specified # of turns and # of particles from a simulation.
        
    Arguments:
        opts (options.Options):  A Synergia options instance containing plotting parameters and path information.
                                 This object also points to the lattice_simulator and other simulation specific details.
        elliptic (Optional[bool]): If True, compute elliptic invariants. Defaults to True.
        
    
    Returns:
        (ndarray): Horizontal stack of arrays of invariant values (H, I)
    
    '''
    
    files = get_file_list(opts)
    twiss = get_sliced_twiss(opts)
    lost = get_lost_particle_list(opts)
    
    if len(lost) > 0:
        #we have lost particles
        opts.lost = lost #store these in opts.lost
        lost = True #make lost a simple flag
    
    hArray = []
    iArray = []
    
    for outfile in files:
        if lost:
            header, particles, lost_particles = get_particles(outfile, True,opts.lost)
        else:
            header, particles = get_particles(outfile, False)
        
        #grab invariants for that file
        hVals, iVals = calc_bunch_H(particles, opts, elliptic)
        
        if hArray == []:   #handle base case
            hArray = np.asarray(hVals)
            iArray = np.asarray(iVals)
        else:
            hArray = np.vstack((hArray,np.asarray(hVals)))
            iArray = np.vstack((iArray,np.asarray(iVals)))
            
    return np.asarray(hArray), np.asarray(iArray)
    
    

def toy_calc_elliptic_Invariant(opts, elliptic=True):
    '''
    Returns both invariants for NL elliptic potential over a specified # of turns and # of particles from a simulation.
     
    The difference with "calc_elliptic_invariant" is that this version uses a fixed Twiss from get_toy_twiss().
       
    Arguments:
        opts (options.Options):  A Synergia options instance containing plotting parameters and path information.
                                 This object also points to the lattice_simulator and other simulation specific details.
    
    Returns:
        hArray, iArray: arrays of invariant values (H, I)
    
    '''
        
    files = get_file_list(opts)
    twiss = get_toy_twiss(opts)
    lost = get_lost_particle_list(opts)
    
    if len(lost) > 0:
        #we have lost particles
        opts.lost = lost #store these in opts.lost
        lost = True #make lost a simple flag
    
    hArray = []
    iArray = []
    
    for outfile in files:
        if lost:
            header, particles, lost_particles = get_particles(outfile, True,opts.lost)
        else:
            header, particles = get_particles(outfile, False)
        
        #grab invariants for that file
        hVals, iVals = calc_bunch_H(particles, opts, elliptic)
        
        if hArray == []:   #handle base case
            hArray = np.asarray(hVals)
            iArray = np.asarray(iVals)
        else:
            hArray = np.vstack((hArray,np.asarray(hVals)))
            iArray = np.vstack((iArray,np.asarray(iVals)))
            
    return np.asarray(hArray), np.asarray(iArray)
    
def toy_calc_Invariant(opts):
    
    '''
    A wrapped for calling toy_calc_elliptic_Invariant with elliptic=False
    
    Arguments:
        opts (options.Options):  A Synergia options instance containing plotting parameters and path information.
                                 This object also points to the lattice_simulator and other simulation specific details.
    
    Returns:
        hArray, iArray: Array of invariant values, and an array of zeroes
    
    '''
    
    return toy_calc_elliptic_Invariant(opts, elliptic=False)


def toy_analyze_invariant_bunch(opts):
    '''
    Computes statistical quantities -- mean, std, variance -- for the invariant data for a given simulations.
    
    This function is designed to further separate particle files containing bunches with different emittances, 
    such that there are opts.macro_particles in each sub bunch having an emittance defined by opts.emits.
    Again, this function relies on get_toy_twiss() and thus cannot be used generically.
    
    Arguments:
        opts (options.Options):  A Synergia options instance containing plotting parameters and path information.
                                 This object also points to the lattice_simulator and other simulation specific details.
    
    Returns:
        hStats, iStats: Arrays of H and I values and relevant statistics, respectively
    
    '''
    
    hArray, iArray = toy_calc_elliptic_Invariant(opts)
    hAfter = hArray[1::]
    iAfter = iArray[1::]
    
    hStats = {}
    iStats = {}
    
    #we assume that the total # of particles is len(opts.emits)*opts.macro_particles
    for index,emit in enumerate(opts.emits):
        
        hVals = []
        iVals = []
        
        #sort through each index here
        for ind in range(index*opts.macro_particles,(index+1)*opts.macro_particles):
            current = hAfter[:,ind]
            
            h_mean = np.mean(current)
            h_std = np.std(current)
            h_variance = (h_std/h_mean)*100.
            hVals.append([ind,h_mean, h_std, h_variance])

        
            #do the same for i
            current = iAfter[:,ind]
            i_mean = np.mean(current)
            i_std = np.std(current)
            i_variance = (i_std/i_mean)*100.
            iVals.append([ind,i_mean, i_std, i_variance])       
        
        
        hStats[str(emit)] = np.asarray(hVals)
        iStats[str(emit)] = np.asarray(iVals)
        
    return hStats, iStats

   
def toy_analyze_invariant(opts):
    '''
    Computes statistical quantities -- mean, std, variance -- for the invariant data for a given simulations.

    Again, this function relies on get_toy_twiss() and thus cannot be used generically.
    
    Arguments:
        opts (options.Options):  A Synergia options instance containing plotting parameters and path information.
                                 This object also points to the lattice_simulator and other simulation specific details.
    
    Returns:
        HS_t, IS_t: Arrays of H and I values and relevant statistics, respectively
    
    '''    
    
    
    hArray, iArray = toy_calc_elliptic_Invariant(opts)
    hAfter = hArray[1::]
    iAfter = iArray[1::]
    
    hStats = []
    iStats = []

    pRange = len(hAfter[0])  #define this particle range to correspond to the total number of particles left
    
    if not pRange == opts.num_total_particles:
           print "Adjusting analysis to account for lost particles"
    
    for ind in range(pRange):
        
        current = hAfter[:,ind] #current returns an array of h values for a given particle # over all turns 
        h_mean = np.mean(current)
        h_std = np.std(current)
        h_variance = (h_std/h_mean)*100.
        
        hStats.append([ind,h_mean, h_std, h_variance])
        
        #do the same for i
        current = iAfter[:,ind]
        i_mean = np.mean(current)
        i_std = np.std(current)
        i_variance = (i_std/i_mean)*100.
 
        iStats.append([ind,i_mean, i_std, i_variance])       
        
    
    #reshape
    HS_l = np.asarray(list(itertools.chain(*hStats)))
    HS_t = np.reshape(HS_l,(-1,4))

    #reshape
    IS_l = np.asarray(list(itertools.chain(*iStats)))
    IS_t = np.reshape(IS_l,(-1,4))
        
    return HS_t, IS_t

def analyze_invariant(opts):
    ''' 
    Computes statistical quantities -- mean, std, variance -- for the invariant data for a given simulations.

    This function grabs the sliced twiss and can be used generically.
    
    Arguments:
        opts (options.Options):  A Synergia options instance containing plotting parameters and path information.
                                 This object also points to the lattice_simulator and other simulation specific details.
    
    Returns:
        HS_t, IS_t: Arrays of H and I values and relevant statistics, respectively
    
    '''
    
    hArray, iArray = calc_elliptic_Invariant(opts)
    hAfter = hArray[1::]
    iAfter = iArray[1::]
    
    hStats = []
    iStats = []

    pRange = len(hAfter[0])  #define this particle range to correspond to the total number of particles left
    
    if not pRange == opts.num_total_particles:
           print "Adjusting analysis to account for lost particles"
    
    for ind in range(pRange):
        
        current = hAfter[:,ind] #current returns an array of h values for a given particle # over all turns 
        h_mean = np.mean(current)
        h_std = np.std(current)
        h_variance = (h_std/h_mean)*100.
        
        hStats.append([ind,h_mean, h_std, h_variance])

        #do the same for i
        current = iAfter[:,ind]
        i_mean = np.mean(current)
        i_std = np.std(current)
        i_variance = (i_std/i_mean)*100.
 
        iStats.append([ind,i_mean, i_std, i_variance])       
        
    
    #reshape
    HS_l = np.asarray(list(itertools.chain(*hStats)))
    HS_t = np.reshape(HS_l,(-1,4))

    #reshape
    IS_l = np.asarray(list(itertools.chain(*iStats)))
    IS_t = np.reshape(IS_l,(-1,4))
        
    return HS_t, IS_t


   
def return_coords(opts):
    '''
    Return the human (default) particle coordinates given an options object pointing to a list of files.
    
    This function grabs the sliced twiss and can be used generically. Note that if opts.plot_lost is specified,
    then one can return the lost particle coordinates instead.
    
    Arguments:
        opts (options.Options):  A Synergia options instance containing plotting parameters and path information.
                                 This object also points to the lattice_simulator and other simulation specific details.
    Returns:
        pcoords (ndarray): Array of coordinates.
    
    '''

    files = get_file_list(opts)
    twiss = get_sliced_twiss(opts)
    lost = get_lost_particle_list(opts)
    
    if not len(lost) == 0:
        #we have lost particles
        opts.lost = lost #store these in opts.lost
        lost = True #make lost a simple flag
    else:
        lost = False
        opts.lost = None
    
    if opts.plot_lost:
         pCoords = get_human_coords(files, True, opts.lost) 
    else:
        pCoords = get_human_coords(files)
    
    return pCoords