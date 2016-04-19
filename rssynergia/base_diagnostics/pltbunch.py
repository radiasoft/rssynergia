import sys
import synergia
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import gridspec

#RC stuff
#mpl.rc('axes', labelsize=12, titlesize=16)


def check_grid_spec():
    '''Constructs a plot illustrating Matplotlib grid spec alignment.
    
    This is a utility function.
    
    '''
    
    fig = plt.figure(figsize=(16,8))
    gs = gridspec.GridSpec(2, 2, width_ratios=[1,1,1])

    ax0 = plt.subplot(gs[0])
    ax0.set_title('gs[0]')

    ax1 = plt.subplot(gs[1])
    ax1.set_title('gs[1]')

    ax2 = plt.subplot(gs[2])
    ax2.set_title('gs[2]')

    ax3 = plt.subplot(gs[3])
    ax3.set_title('gs[3]')

    plt.show()


def plot_bunch(bunch):
    '''
    Plot the coordinate space and x, y phase spaces of a Synergia bunch.
    
    List of plots: x-y, x-x', y-y'
    
    Arguments:
        bunch (synergia.bunch.bunch): A Synergia bunch object.
        
    When used with an IPython backend, figure should be displayed.
        
    '''
    
    #get particle properties from bunch
    myParticles = bunch.get_local_particles()
    xp = myParticles[:,bunch.xp]
    x = myParticles[:,bunch.x]
    yp = myParticles[:,bunch.yp]
    y = myParticles[:,bunch.y]
    
    xmax = np.max(x)
    xpmax = np.max(xp)
    ymax = np.max(y)
    ypmax = np.max(yp)
    #one way to use subplots
    #fig, (ax0, ax1, ax2)  = plt.subplots(ncols=3, figsize=(10,6))

    #another way - use gridspec
    fig = plt.figure(figsize=(15,5))
    gs = gridspec.GridSpec(1, 3, width_ratios=[1,1,1]) 
    
    ax0 = plt.subplot(gs[0])
    ax0.scatter(x, y, c='k',s=16)
    ax0.set_title('X-Y Coordinate Space',fontsize='16')
    ax0.set_xlim([-1.5*xmax,1.5*xmax])
    ax0.set_ylim([-1.5*ymax,1.5*ymax])
    #ax0.set_aspect(aspect=2.0)
    
    ax1 = plt.subplot(gs[1])
    ax1.scatter(x, xp, c='b',s=16)
    ax1.set_title('X Trace Space',fontsize='16')
    ax1.set_xlim([-1.5*xmax,1.5*xmax])
    ax1.set_ylim([-1.5*xpmax,1.5*xpmax])
    #ax1.set_aspect(aspect=2.0)
    
    ax2 = plt.subplot(gs[2])
    ax2.scatter(y, yp, c='r',s=16)
    ax2.set_title('Y Trace Space',fontsize='16')
    ax2.set_xlim([-1.5*ymax,1.5*ymax])
    ax2.set_ylim([-1.5*ypmax,1.5*ypmax])
    #ax2.set_aspect(aspect=2.0)
    
    # Tweak spacing between subplots to prevent labels from overlapping
    #plt.subplots_adjust(hspace=2)
    
    #set figure title
    fig.canvas.set_window_title('Synergia Phase Space Distribution')
    fig.tight_layout()
    plt.show()
    

def plot_long(bunch):
    '''
    Plot the longitudinal coordinate and phase spaces of a Synergia bunch.

    List of plots: x-y, z-z', z-y
    
    Arguments:
        bunch (synergia.bunch.bunch): A Synergia bunch object.
        
    When used with an IPython backend, figure should be displayed.
    '''
    
    #get particle properties from bunch
    myParticles = bunch.get_local_particles()
    x = myParticles[:,bunch.x]
    y = myParticles[:,bunch.y]
    z = myParticles[:,bunch.z]
    zp = myParticles[:,bunch.zp]

    xmax = np.max(x)
    ymax = np.max(y)
    zmax = np.max(z)
    zpmax = np.max(zp)

    
    #another way - use gridspec
    fig = plt.figure(figsize=(15,5))
    gs = gridspec.GridSpec(1, 3, width_ratios=[1,1,1]) 
    
    ax0 = plt.subplot(gs[0])
    ax0.scatter(x, y, c='k',s=16)
    ax0.set_title('X-Y Coordinate Space',fontsize='16')
    ax0.set_xlim([-1.5*xmax,1.5*xmax])
    ax0.set_ylim([-1.5*ymax,1.5*ymax])
    #ax0.set_aspect(aspect=2.0)
    
    ax1 = plt.subplot(gs[1])
    ax1.scatter(z, zp, c='b',s=16)
    ax1.set_title('Z Trace Space',fontsize='16')
    ax1.set_xlim([-1.5*zmax,1.5*zmax])
    ax1.set_ylim([-1.5*zpmax,1.5*zpmax])
    #ax1.set_aspect(aspect=2.0)
    
    ax2 = plt.subplot(gs[2])
    ax2.scatter(z, y, c='r',s=16)
    ax2.set_title('Z-Y Coordinate Space',fontsize='16')
    ax2.set_xlim([-1.5*zmax,1.5*zmax])
    ax2.set_ylim([-1.5*ymax,1.5*ymax])


    #ax2.set_aspect(aspect=2.0)
    
    # Tweak spacing between subplots to prevent labels from overlapping
    #plt.subplots_adjust(hspace=2)
    
    #set figure title
    fig.canvas.set_window_title('Synergia Phase Space Distribution')
    fig.tight_layout()
    plt.show()



def plt_bunch_dpop(bunch,opts,lost=None):
    '''A plotting script displaying different phase space plots - able to be split up by dpop
    
    List of plots: x-y, x'-y', x-x', y-y'
    
    Arguments:
        bunch (synergia.bunch.bunch): A Synergia bunch object.
        opts (options.Options): A Synergia options instance specifying dpop list, lattice name, turn #, etc.
        lost (Optional[bool]): If 
        
    When used with an IPython backend, figure should be displayed. If opts.save is specified,
    then it should write each plot to a file.    
    
    '''
    
    
    if lost:
        full_particles = bunch.get_local_particles()
        myParticles = full_particles[lost] #extract only the ones with indices from lost
    else:
        myParticles = bunch.get_local_particles()
    
    
    coord_max = 2.*max([np.max(np.abs(myParticles[:,0])),np.abs(np.max(myParticles[:,2]))])
    mom_max  = 2.*max([np.max(np.abs(myParticles[:,1])),np.max(np.abs(myParticles[:,3]))])
    
    #for index,val in enumerate(opts.dpops):
    for index,val in enumerate(opts.dpops):
        st = index*opts.n_macro
        en = (index+1)*opts.n_macro - 1

        particles = myParticles[st:en]
        
        xp = particles[:,bunch.xp]
        x = particles[:,bunch.x]
        yp = particles[:,bunch.yp]
        y = particles[:,bunch.y]
        
        #use gridspec
        fig = plt.figure(figsize=(12,8))
        gs = gridspec.GridSpec(2, 2, width_ratios=[1,1,1]) 
        
        ax0 = plt.subplot(gs[0])
        ax0.scatter(x, y, c='b', s=1)
        ax0.set_title('x-y coordinate space')
        ax0.set_xlim([-1*coord_max,coord_max])
        ax0.set_ylim([-1*coord_max,coord_max])
    
        ax1 = plt.subplot(gs[1])
        ax1.scatter(xp, yp, c='r', s=1)
        ax1.set_title('px - py momentum space')
        ax1.set_xlim([-1*mom_max,mom_max])
        ax1.set_ylim([-1*mom_max,mom_max])
    
        ax2 = plt.subplot(gs[2])
        ax2.scatter(x, xp, c='g', s=1)
        ax2.set_title('x trace space')
        ax2.set_xlim([-1*coord_max,coord_max])
        ax2.set_ylim([-1*mom_max,mom_max])
           
        ax3 = plt.subplot(gs[3])
        ax3.scatter(y, yp, c='k', s=1)
        ax3.set_title('y trace space')
        ax3.set_xlim([-1*coord_max,coord_max])
        ax3.set_ylim([-1*mom_max,mom_max])       
        plt.show()    

        if opts.save:
            #save each dpop slice
            sv_title = 'Phase_Space'+'_' + str(val)+'_'+ str(opts.turns) + '_turns_'+  opts.lattice_name + '.pdf'
            fig.savefig(sv_title, bbox_inches='tight') 

    