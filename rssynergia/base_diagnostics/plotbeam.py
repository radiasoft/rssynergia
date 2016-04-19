#!/usr/bin/env python

import sys
import tables
import numpy as np
import math
import matplotlib.pyplot as plt
import options

coords = {}
coords['x'] = 0
coords['xp'] = 1
coords['y'] = 2
coords['yp'] = 3
coords['cdt'] = 4
coords['dpop'] = 5
coords['id'] = 6


def plotbeam(opts, particles, header):

    '''Returns a 2D plot of particle distribution in the choosen coordinate space
    
    Arguments:
        opts (options.Options): A Synergia options instance 
        particles (ndarray): An array of particles, organized according to the coords dictionary
        header (dict): A Python dictionary with metadata for the particle array
    
    '''

    h = particles[:,coords[opts.hcoord]]
    v = particles[:,coords[opts.vcoord]]
    #vcoords = particles[:,coords[vc]] for vc in opts.vcoord

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
    
    '''
    Reads an input file and returns a numpy array of particles and a dictionary of root values
    
    Arguments:
        opts (options.Options): A Synergia options instance
    
    '''
    
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
    '''
    Plot a beam of particles given an options object input.
    
    Arguments:
        opts (options.Options): A Synergia options instance
    
    '''
    
    header, particles = get_particles(opts)

    if opts.plots ==1:
        plotbeam(opts, particles,header)
        