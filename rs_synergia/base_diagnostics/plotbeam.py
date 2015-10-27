#!/usr/bin/env python

import sys
import tables
import numpy as np
import math
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


def plot_one(opts, particles,header):

    h = particles[:,coords[opts.hcoord]]
    v = particles[:,coords[opts.vcoord]]

    fig = plt.figure()
    ax = plt.gca()
    
    ax.scatter(h,v, c ='b')
    ax.set_aspect('equal', 'datalim')
    #plt.plot(h,v, 'o')
    plt.xlabel(opts.hcoord,fontsize=12)
    plt.ylabel(opts.vcoord,fontsize=12)
    title = 'Particle distribution at s = ' + str(header['t_len'])
    if not opts.lattice_name== None:
        title = title + ' for lattice ' + opts.lattice_name
    plt.title(title, y=1.05, fontsize=14)
    plt.show()
    
    if opts.save:
        sv_title = 'Beam_' + opts.lattice_name + '.pdf'
        fig.savefig(sv_title, bbox_inches='tight')
    
    
#define an option to replicate the pltbunch.plot_bunch function?


def get_particles(opts):
    
    '''Reads an input file and returns a numpy array of particles and a dictionary of root values'''
    
    f = tables.openFile(opts.inputfile, 'r')
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
    
    return header,particles
    

def plot_beam(opts):
    
    
    header, particles = get_particles(opts)

    if opts.plots ==1:
        plot_one(opts, particles,header)
        