#A collection of miscellaneous useful scripts. These will be refined and reorganized in time.
#Author: Nathan Cook
#10/10/2015

import numpy as np
import scipy.constants as const
import synergia
import os
import matplotlib.pyplot as plt
import tables
from matplotlib import gridspec

def params(mass,kE):
    '''
    Return dictionary of values relating momentum, gamma, etc for a given input mass and energy
    
    Args:
        mass (float): mass of the particle in GeV
        kE (float): kinetic energy of particle
    
    
    '''
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
    '''Returns the rms beam size in a given dimension
    
    Args:
        dim (str): dimension of calculation ('x','y')
        particles (ndarray): array of particles of the form (x,px,y,py)
    
    Returns:
        sig_w (float): RMS value of beam envelop in the desired dimension
    
    '''
    
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
        
    return  sig_w 
    
    
def get_emittance(dim, bunch):
    '''Returns the statistical emittance (in x or y) of a bunch  
    centered at (0,0) and not rotated with respect to the x-y axes.

    Args:
        dim (str): dimension of calculation ('x','y')
        particles (ndarray): array of particles of the form (x,px,y,py)
    
    Returns:
        emittance (float): RMS emittance value in the desired dimension   
    
    '''
    
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
    
    return  sig_w * sig_w_prime


def get_normalized_emittance(dim, bunch, beta, gamma):
    '''Return the normalized emittance of a bunch given relativistic beta and gamma
    
    Args:
        dim (str): dimension of calculation ('x','y')
        bunch (ndarray): array of particles of the form (x,px,y,py)
        beta (float): rleativistic beta of bunch
        gamma (float): relativistic gamma of bunch
    
    Returns:
        n_emit (float): RMS normalized emittance value in the desired dimension       

    '''
    
    return get_emittance(dim, bunch)*(beta*gamma)
    
def calc_geometric_emittance(norm_emit,beta,gamma):
    '''Return the geometric emittance given a normalized emittance value and a beta/gamma.

    Args:
        norm_emit (float): a normalized emittance value
        beta (float): rleativistic beta of bunch
        gamma (float): relativistic gamma of bunch
    
    Returns:
        g_emit (float): RMS geometric emittance value
        
    '''
    
    return norm_emit/(beta*gamma)
    
def calc_normalized_emittance(g_emit,beta,gamma):
    '''Return the normalized emittance given a geometric emittance value and a beta/gamma.
    
    Args:
        g_emit (float): a geometric emittance value
        beta (float): rleativistic beta of bunch
        gamma (float): relativistic gamma of bunch
    
    Returns:
        norm_emit (float): RMS normalized emittance value
    
    '''
    
    return g_emit*(beta*gamma)


def calc_properties(bunch,ref):
    '''Print a list of different bunch properties
    
    Args:
        bunch (synergia.bunch.bunch): a Synergia bunch object
        ref (synergia.foundation.reference_particle): a Synergia reference particle object
    
    '''
    
    beta = ref.get_beta()
    gamma = ref.get_gamma()
    
    g_emit_x = get_emittance('x',bunch)
    g_emit_y = get_emittance('y',bunch)
    particles = bunch.get_local_particles()
    x_vals = particles[:,0]
    y_vals = particles[:,2]
    xp_vals = particles[:,1]
    yp_vals = particles[:,3]
    
    xenv = np.sqrt(np.mean(x_vals**2))
    print "rms envelope x: {} mm".format(xenv*1.e3)
    
    yenv = np.sqrt(np.mean(y_vals**2))
    print "rms envelope y: {} mm".format(yenv*1.e3)
    
    xmax = np.max(x_vals)
    print "maximum x value is : {} mm".format(xmax*1.e3)
    
    ymax = np.max(y_vals)
    print "maximum y value is : {} mm".format(ymax*1.e3)
    
    xbet = np.mean(x_vals**2)/g_emit_x
    print "rms beta x: {}".format(xbet)
    
    ybet = np.mean(y_vals**2)/g_emit_y
    print "rms beta y: {}".format(ybet)
    
    emitx = get_emittance('x',bunch)
    print "geometric emittance x: {} mm-mrad".format(emitx*1.e6)
    
    emity = get_emittance('y',bunch)
    print "geometric emittance y: {} mm-mrad".format(emity*1.e6)    
        
    n_emitx = get_normalized_emittance('x',bunch,beta,gamma)
    print "normalized emittance x: {} mm-mrad".format(n_emitx*1.e6)
    
    n_emity = get_normalized_emittance('y',bunch,beta,gamma)
    print "normalized emittance y: {} mm-mrad".format(n_emity*1.e6)
    
    x2_mean = np.mean(xp_vals**2)
    print "mean of xp^2 : {}".format(x2_mean)
    
    y2_mean = np.mean(yp_vals**2)
    print "mean of yp^2 : {}".format(y2_mean)
    
    emit_total_x = xmax**2/xbet
    print "total geometric emittance x: {} mm-mrad".format(emit_total_x*1.e6)
    
    emit_total_y = ymax**2/ybet
    print "total geometric emittance y: {} mm-mrad".format(emit_total_y*1.e6)


def get_base_nll(l0, mu0, t, c):
    '''Construct the nonlinear element. Taken from madx script by A. Valishev. 
    Verified by David's python script for pyOrbit.
    
    Args:
        l0 (float): length of nonlinear drift section
        mu0 (float): tune advance (phase/2pi) through the drift section
        t (float): nonlinear magnet strength parameter
        c (float): nonlinear magnet aperture parameter
    
    Returns:
        f0, betae, alphae, betas: array of values - focal length, entrance beta, entrance alpha,
                                  beta at middle of NL drift
    
    '''

    f0=l0/4.0*(1.0+1.0/np.tan(np.pi*mu0)**2); #focal length
    betae=l0/np.sqrt(1.0-(1.0-l0/2.0/f0)**2); #entrance beta function
    alphae=l0/2.0/f0/np.sqrt(1.0-(1.0-l0/2.0/f0)**2); #entrance alpha function
    betas=l0*(1-l0/4.0/f0)/np.sqrt(1.0-(1.0-l0/2.0/f0)**2); #middle beta function
    
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
