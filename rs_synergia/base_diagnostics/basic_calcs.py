#A collection of miscellaneous useful scripts. These will be refined and reorganized in time.

import numpy as np
import scipy.constants as const
import synergia
import os
import matplotlib.pyplot as plt

def params(mass,kE):
    '''return dictionary of values relating momentum, gamma, etc for a given input mass and energy'''
    vals = {}
    vals['units'] = 'GeV'
    vals['mass'] = mass
    vals['kE'] = kE
    
    gamma = np.divide(kE,mass) + 1
    total_energy = gamma*mass #total energy in GeV
    vals['gamma'] = gamma
    vals['total_energy'] = total_energy
    vals['beta'] = np.sqrt(1. - 1./gamma**2)
    vals['pc'] = np.sqrt(total_energy**2 - mass**2)
    
    return vals
    
    
def get_rms_envelope(dim, particles):
    '''Returns the rms beam size in a given dimension'''
    
    if dim == 'x':
        ind = [0,1]
    elif dim == 'y':
        ind = [2,3]
    else:
        print "Please specify 'x' or 'y' axis"
  
    #get particles and assign w, w' 
    w = particles[:,ind[0]]
    p_num = w.shape[0]
    sig_w = np.sqrt(np.sum(w**2)/p_num)
    
    #w_prime = particles[:,ind[1]]
    #sig_w_prime = np.sqrt(np.sum(w_prime**2)/p_num)
    
    #there could be a factor of 2 here that needs to be included, depending on convention!
    
    return  sig_w 
    


    
def get_emittance(dim, bunch):
    '''Returns the statistical emittance (in x or y) of a bunch which centered at (0,0) and not rotated with respect to the x-y axes'''
    
    if dim == 'x':
        ind = [0,1]
    elif dim == 'y':
        ind = [2,3]
    else:
        print "Please specify 'x' or 'y' axis"
  
    #get particles and assign w, w' 
    particles = bunch.get_local_particles()
    w = particles[:,ind[0]]
    w_prime = particles[:,ind[1]]
    
    p_num = w.shape[0]
    sig_w = np.sqrt(np.sum(w**2)/p_num)
    
    sig_w_prime = np.sqrt(np.sum(w_prime**2)/p_num)
    
    #there could be a factor of 2 here that needs to be included, depending on convention!
    
    return  sig_w * sig_w_prime


def get_normalized_emittance(dim, bunch, beta, gamma):
    '''Calculate the normalized emittance for a bunch with a given gamma and beta'''
    
    return get_emittance(dim, bunch)*(beta*gamma)
    
def get_geometric_emittance(norm_emit,beta,gamma):
    '''Return the geometric emittance given a normalized emittance value and a beta/gamma.'''
    
    return norm_emit/(beta*gamma)
    
    
def get_base_nll(nn, l0, mu0, t, c):
    '''Construct the nonlinear element. Taken from madx script by A. Valishev. 
    Verified by David's python script for pyOrbit.
    
    '''
    #musect=mu0+0.5;
    f0=l0/4.0*(1.0+1.0/np.tan(np.pi*mu0)**2); #focal length
    betae=l0/np.sqrt(1.0-(1.0-l0/2.0/f0)**2); #entrance beta function
    alphae=l0/2.0/f0/np.sqrt(1.0-(1.0-l0/2.0/f0)**2); #entrance alpha function
    betas=l0*(1-l0/4.0/f0)/np.sqrt(1.0-(1.0-l0/2.0/f0)**2); #middle beta function
    #value,f0,betae,alfae,betas;
    
    return [f0, betae,alphae,betas]


def make_Hvals(opts, ID=0):
    
    Hvals = {}
    Ivals = {}
    
    HArray, IArray = elliptic_sp.toy_calc_Invariant(opts)
    
    for index,emit in enumerate(opts.emits):
        
        ind = int(opts.macro_particles*index + ID)
        Hvals[str(emit)] = HArray[1:,ind]*1.e6
        Ivals[str(emit)] = IArray[1:,ind]*1.e6
        
    return Hvals, Ivals


import tables
from matplotlib import gridspec

def plot_bunch(particles_file, opts=None, part=None):
    '''Returns a synergia bunch object constructed from a particles output file.'''
    
    if opts:
        inputfile = os.path.join(opts.output_dir,particles_file)
    else:
        inputfile = particles_file
    f = tables.openFile(inputfile, 'r')
    particles = f.root.particles.read()

    npart = particles.shape[0]
    
    if part:  
        xp = particles[part[0]:part[1],1]
        x = particles[part[0]:part[1],0]
        yp = particles[part[0]:part[1],3]
        y = particles[part[0]:part[1],2]
    else:
        xp = particles[:,1]
        x = particles[:,0]
        yp = particles[:,3]
        y = particles[:,2]        
    
    #one way to use subplots
    #fig, (ax0, ax1, ax2)  = plt.subplots(ncols=3, figsize=(10,6))

    #another way - use gridspec
    #fig = plt.figure(figsize=(9.9,3.3))
    scaleFig=4.5
    fig = plt.figure(figsize=(3*scaleFig,1*scaleFig))
    gs = gridspec.GridSpec(1, 3, width_ratios=[1,1,1]) 
    
    ax0 = plt.subplot(gs[0])
    ax0.scatter(x, y, c='b')
    ax0.set_title('X-Y Coordinate Space',fontsize=14)
    ax0.set_xlabel('x [m]')
    ax0.set_ylabel('y [m]')
    #ax0.set_aspect(aspect=2.0)
    
    ax1 = plt.subplot(gs[1])
    ax1.scatter(x, xp, c='r')
    ax1.set_title('X Trace Space',fontsize=14)
    ax1.set_xlabel('x [m]')
    ax1.set_ylabel('x-prime [rad]')
    #ax1.set_aspect(aspect=2.0)
    
    ax2 = plt.subplot(gs[2])
    ax2.scatter(y, yp, c='g')
    ax2.set_title('Y Trace Space',fontsize=14)
    ax2.set_xlabel('y [m]')
    ax2.set_ylabel('y-prime [rad]')
    #ax2.set_aspect(aspect=2.0)
    
    
    ax2.tick_params(axis='x', pad=5)
    
    # Tweak spacing between subplots to prevent labels from overlapping
    #plt.subplots_adjust(hspace=2)
    
    #set figure title
    #fig.canvas.set_window_title('Bunch Distribution')
    fig.tight_layout()
    
    #tn = f.root.tlen[()] #cumulative distance travelled
    sn = f.root.s_n[()] #s value along lattice for a given turn (e.g. tn mod lattice_length)
    turn_num = f.root.rep[()]
    
    if not sn == 0:
        turn_num +=1
    
    svtitle = "bunch_turn_{}".format(turn_num)
    supertitle = 'Bunch distribution for turn {}'.format(turn_num+1)
    
    if opts:
        svtitle = svtitle +'_' + opts.name + '.pdf'
        supertitle = 'Bunch distribution for NL tune {} at turn {}'.format(opts.new_tune,turn_num)
        #supertitle += ' with NL tune {}'.format(opts.new_tune)
    else:
        svtitle += '.pdf'
    
    plt.suptitle(supertitle, y=1+0.07, fontsize=18)
    plt.show()
    
    fig.savefig(svtitle, bbox_inches='tight')
    print "Saving plot to",svtitle
