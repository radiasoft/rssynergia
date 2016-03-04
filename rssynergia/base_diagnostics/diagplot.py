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
    def __init__(self, label, x_attr, y_attr, y_index1=None, y_index2=None):
        self.label = label
        self.x_attr = x_attr
        self.y_attr = y_attr
        self.y_index1 = y_index1
        self.y_index2 = y_index2
        

#generate_plotparams constructs the associated values capable of being plotted
def generate_plotparams():
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


def plot_turn(x,yVals,start,names, numE, lname, emit, save, betas,length=0):
    '''Makes plot of x vs. yVals and saves it.'''
    
    #if the length is specified, use it as the intervl descriptor
    if not (length==0):
        #plot the entire thing since length parameter truncates plot
        st = 0
        et = len(x)
    else:
        #just plot the specified interval  
        interval = numE # number of steps per turn

        st = start*interval
        et = (start+1)*interval-1

    #handle multiple y arrays
    for index,(name,y) in enumerate(map(None,names,yVals)):
        
        #plt.plot(x[st:et],y[st:et], '-', label=name[1:])
        
        #newy is the plot holder
        newy = y   
        #deal with betas
        if betas and name == '_x_std':
            newy = makeBeta(y,emit[0])
            names[index] = '_beta_x' #add underscore b/c the label will be cut
            #yVals[index] = newy
        elif betas and name == '_y_std':
            newy = makeBeta(y,emit[1])
            names[index] = '_beta_y'
            #y = newy
            #yVals[index] = newy
            
        #print names       
        
        #add each plot with correct label
        plt.plot(x[st:et],newy[st:et], '-', label=names[index][1:])
    
    #print names
        
    fig1 = plt.gcf() #get current figure handle for saving purposes
    plt.legend(loc='best')
    plt.xlabel('s', fontsize=12)
    
    yNames = ', '.join([n[1:] for n in names])
    plt.ylabel(yNames, fontsize=12)
    
    #limit x-axis according to 'length' specification
    if not (length==0):
        axes = plt.gca()
        axes.set_xlim([0,length])
        if betas:
            name = lname+''.join(names)+'_betas.pdf'
            title_name = 'Estimated RMS betas for lattice ' + lname  
        else:
            name = lname+''.join(names)+'.pdf'
            title_name = yNames + ' for lattice ' + lname
    
    
    #plotting based on turn
    else:
        if betas:
            #names is a list so have to join with ''
            name = lname+''.join(names)+'_turn'+str(start+1)+'_betas.pdf'
            title_name = 'Estimated RMS betas for lattice ' + lname+': Turn '+str(start+1)
        
        else:
            name = lname+''.join(names)+'_turn'+str(start+1)+'.pdf'
            title_name = yNames + ' for lattice' + lname+': Turn'+str(start+1)
    
    plt.title(title_name, y=1.05, fontsize=15)    
    plt.show()    
    if save:
        fig1.savefig(name)


def makeBeta(vals, emit):
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

def diagPlot(filename, plots, numE,lname, start=0, save=False, betas=False, emit=[1.0e-6,1.0e-6], length=0):
    '''Runs the plotting function to create a plot
    
    Arguments:
        -filename - '.h5' file at minimum basic diagnostics
        -plots - array of y-Axis values to plot
        -numE - number of steps per turn
        -lname - lattice name
        -start - turn start
        -emit - a vector of emittance values - defaults to 1e10-6
        -betas (optional) - will roughly estimate corresponding betas for linear lattice
        -length (optional) - will fix the lattice length based on this # rather than the lattice length
        -save (optional) - will save the resulting plot  - defaults to false
    
    '''

    #get plot values dictionary
    plotVals = getPlotVals(filename, plots)
    
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
    
    plot_turn(xmaster,yVals,start,nms, numE, lname, emit, save, betas, length)