#basic_diag_plot.py
import matplotlib.pyplot as plt
import numpy as np
import tables
import synergia
import sys

#define coordinates dictionary
coords = {}
coords['x'] = 0
coords['xp'] = 1
coords['y'] = 2
coords['yp'] = 3
coords['z'] = 4
coords['zp'] = 5

#Define Params class
class Params:
    '''Basic class for coordinating attributes for Synergia data for plotting and manipulation.
    Useful as a lookup for pairing diagnostics output.
    
    Args:
        label (str): a label for the parameter in consideration
        x_attr (str): the x-axis value for plotting (e.g. s)
        y_attr (str): the y- axis value for plotting (e.g. emitx)
        y_index1 (Optional[int]): index of y_attr in particle or diagnostics array. Defaults to None.
        y_index2 (Optional[int]): second index of y_attr in particle or diagnostics array. Defaults to None. 
    '''
    
    def __init__(self, label, x_attr, y_attr, y_index1=None, y_index2=None):
        self.label = label
        self.x_attr = x_attr
        self.y_attr = y_attr
        self.y_index1 = y_index1
        self.y_index2 = y_index2
        

#generate_plotparams constructs the associated values capable of being plotted
def generate_plotparams():
    '''
    Constructs a fixed dictionary of plotparams objects coordinating plottable values from diagnostics.
    
    '''
    
    
    plotparams = {}
    for label in coords.keys():
        for label2 in coords.keys():
            if coords[label2] > coords[label]:
                corr = label + '_' + label2 + '_corr'
                plotparams[corr] = Params(corr, 's', 'corr',
                                           coords[label], coords[label2])
                mom2 = label + '_' + label2 + '_mom2'
                plotparams[mom2] = Params(mom2, 's', 'mom2',
                                           coords[label], coords[label2])
        std = label + '_std'
        plotparams[std] = Params(std, 's', 'std', coords[label])
        mean = label + '_mean'
        plotparams[mean] = Params(mean, 's', 'mean', coords[label])
    plotparams['x_emit'] = Params('x_emit', 's', 'emitx')
    plotparams['y_emit'] = Params('y_emit', 's', 'emity')
    plotparams['z_emit'] = Params('z_emit', 's', 'emitz')
    plotparams['xy_emit'] = Params('xy_emit', 's', 'emitxy')
    plotparams['xyz_emit'] = Params('xyz_emit', 's', 'emitxyz')
    plotparams['particles'] = Params('particles', 's', 'num_particles')
    return plotparams


def plot_turn(x,yVals,opts,names,betas,length=None):
    '''Makes plot of x vs. yVals and saves it.
    
    Args:
        x (array-like): array of x-axis values
        y (array-like): array of y-axis values
        opts (options instance): Synergia options instance with metadata for simulation
        names (list):  a list of names of the coordinates is being plotted
        length (float): Sets a maximum on the x-axis
    
    '''
    
    #if the length is specified, use it as the interval descriptor
    if length:
        #plot the entire thing since length parameter truncates plot
        st = 0
        et = len(x)
    else:
        #just plot the specified interval  
        #interval = len(opts.lattice_simulator.get_slices()) # number of steps per turn
        interval = opts.steps
        interval = 896
        st = opts.start*interval
        
        if opts.turns:
            et = (opts.start+opts.turns)*interval-1
        else:
            et = (opts.start+1)*interval-1

    fig1 = plt.figure(figsize=(12,8))
    
    #fig1 = plt.gcf() #get current figure handle for saving purposes
    plt.xlabel('s [m]', fontsize=14)
    
    yNames = ', '.join([n[1:] for n in names])
    names2 = []
    
    for n in names:
        if n.endswith('std'):
            names2.append('$\sigma_{}$'.format(n[1]))

    yNames = ', '.join([n for n in names2])
    plt.ylabel(yNames + ' [mm]', fontsize=14)
    
    #handle multiple y arrays
    for index,(name,y) in enumerate(map(None,names,yVals)):
        
        #plt.plot(x[st:et],y[st:et], '-', label=name[1:])
        
        #newy is the plot holder
        newy = y   

        #add each plot with correct label
        plt.plot(x-opts.lattice.get_length()*(opts.start),newy*1000, '-', label=names2[index])

    ax = plt.gca()
    
    if not opts.turns:
        ax.set_xlim([0,opts.lattice.get_length()])
    else:
       ax.set_xlim([0,opts.turns*opts.lattice.get_length()])

    if length:
        axes = plt.gca()
        axes.set_xlim([0,length])
        name = opts.lattice_name+''.join(names)+'.pdf'
        title_name = yNames + ' for lattice ' + opts.lattice_name
    
    #plotting based on turn
    elif opts.turns:
        sv_name = opts.lattice_name+''.join(names)+'_turn'+str(opts.start)+'-'+str(opts.start+opts.turns)+'.pdf'
        title_name = yNames + ' for lattice ' + opts.lattice_name+': Turns '+str(opts.start)+'-'+str(opts.start+opts.turns)
            
    #just plot 1 turn based on the starting position
    else:
        sv_name = opts.lattice_name+''.join(names)+'_turn'+str(opts.start+1)+'.pdf'
        title_name = yNames + ' for lattice ' + opts.lattice_name+': Turn '+str(opts.start+1)
    
    plt.legend(loc='best', fontsize=16)
    plt.title(title_name, y=1.02, fontsize=16)    
    plt.show()    
    if opts.save:
        fig1.savefig(sv_name, bbox_inches='tight')


def makeBeta(vals, emit):
    '''Returns a quick approximation of beta functions given an emittance and assuming a linear lattice'''
    #for RMS, sigma is sqrt(beta*epsilon), so beta is sigma^2/epsilon
    return np.multiply(vals,vals)/emit

def getPlotVals(filename, plots):

    #plots is a list of names of variables for plotting
    plotVals = {}

    #generate plot params object
    plotparams = generate_plotparams()
    
    f = tables.openFile(filename, 'r')
    
    #assume all x values are master 's' for now
    param1 = plotparams[plots[0]]
    x = getattr(f.root, param1.x_attr).read()
    plotVals[param1.x_attr] = x
    
    for plot in plots:
    
        params = plotparams[plot]
        #params.y_attr returns a 6D vector for 1st order moments, and a 6x6 matrix for 2nd order
        ymaster = getattr(f.root, params.y_attr).read()
        #specify just the value desired using the y_index1 attribute
        y = ymaster[params.y_index1]
        #put the arrays in the plotVals dictionary
        plotVals[params.label] = y
    
    f.close()
    
    return plotVals

def diagPlot(opts, start=None, length=None):
    '''Runs the plot_turn function to create a plot
    
    Arguments:
        opts (options.Options): An instance of the Synergia options class -Must specify: input_file, lattice_name, lattice_simulator, start (starting turn #)
        length (Optional[float]): fixes the endpoint of the plot (default is lattice length)
    '''

    #get plot values dictionary
    #plotVals = getPlotVals(filename, plots)
    plotVals = getPlotVals(opts.inputfile, opts.plots)
    
    #construct yVals list
    yVals = []
    #construct names list
    nms = []
    #counting index
    index = 0
    
    for key in plotVals.keys():
        #add name to names list and y-values to yvalues list
        nm = ''
        if not(key == 's'):
            nm += '_'
            nm += key
            nms.append(nm)
            yVals.append(plotVals[key])
            index += 1
        
    #just call for the first x variable - all should be parameterized by 's' for now
    xmaster = plotVals['s']
    
    plot_turn(xmaster,yVals,opts,nms, length)