#poincare.py
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


def getCoordVals(filename, plots):
    '''Return an array [xVals, yVals] for the given plot specifications and filename '''
    
    f = tables.openFile(filename, 'r') 
    ## for tracking coordinates, load via get_attr('coords')
    particle_coords = getattr(f.root, "coords").read()
    #particle_coords contains the 6D coordinate vector at all tracked positions
    f.close()
    
    xVals = particle_coords[coords[plots[0]], :]
    yVals = particle_coords[coords[plots[1]], :]
    
    return [xVals,yVals]



def plotP(xC,yC, start, end, names, steps, lname,save):
    '''Make a poincare plot for tracked particles'''
    
    # number of steps per turn
    #interval = steps_per_element*len(lattice.get_elements())
    interval = steps

    st = start*interval
    et = end*interval-1

    #print str(et-st) + " total steps plotted"

    xName = names[0]
    yName = names[1]
    
    #add each plot with correct label
    p = plt.plot(xC[st:et],yC[st:et], 'o', label=yName)
    #plt.plot(x, y, 'o')
    plt.setp(p, "markeredgewidth", 0)

    plt.legend(loc='best')
    plt.legend(numpoints=1)
    plt.xlabel(xName)
    plt.ylabel(yName)
    
    #if one turn plotted, ajust title
    if end-start == 1:
        name = lname+'-'.join(names)+'_turn'+str(start+1)+'.pdf'
        title_name = '-'.join(names) + ' for ' + lname+': Turn '+str(start+1)
    else:
        name = lname+'-'.join(names)+'_turns'+str(start+1)+'-'+str(end)+'.pdf'
        title_name = '-'.join(names) + ' for ' + lname+': Turns '+str(start+1)+'-'+str(end)       
    
    plt.title(title_name)
    fig1 = plt.gcf() #get current figure handle for plotting purposes
    plt.show()
    if save==True:    
        fig1.savefig(name)

    
def poincarePlot(filename, xAx, yAx, numE, lname, start=0, end=1, save=False):
    '''Create and save a one-turn(or more) poincare plot as desired
    
    Arguments:
        -filename - '.h5' file at minimum basic diagnostics
        -xAx,yAx - coordinates to be plotted
        -numE - number of steps per turn
        -lname - lattice name
        -start - turn start - defaults to 0
        -end - turn end - defaults to 1
    
    '''

    #get particle coordinates
    plotCoords = getCoordVals(filename, [xAx,yAx])
    
    #construct names
    nms = [xAx,yAx]
    
    plotP(plotCoords[0],plotCoords[1],start,end,nms,numE,lname,save)