import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import synergia
import math
import matplotlib
import options



########################################### Plotting Functions ###############################################

def lf_plot(lf_info, opts):
    '''
        Plot lattice functions on same plot
        
        Arguments:
        lf_info - array of lattice element dictionaries, each containing lattice functions
        opts - Options object with plot options (including lattice, lattice_simulator, etc)

    '''

    #general plot settings
    matplotlib.rcParams['figure.autolayout'] = True

    
    #get s coordinate and put in numpy array
    ss = np.array([lfd['s'] for lfd in lf_info])
    
    #define figure
    fig1 = plt.figure(figsize=(12,8))
    
    #grab lattice functions as needed
    fn_labels = [] #keep track of modified function labels
    for fn in opts.lf_fns:
        #create y array for each lattice function
        y = np.array([lfd[fn] for lfd in lf_info])
        #add to plot
        ax = fig1.gca()
        if fn.startswith('beta') or fn.startswith('alpha'):
            fn_label = '$\{}$'.format(fn)
        else:
            fn_label = '${}$'.format(fn)
        fn_labels.append(fn_label)
        ax.plot(ss, y, label=fn_label)
        
    
    #add plot labels
    plt.xlabel('s', fontsize=18)
    plt.ylabel(', '.join(fn_labels), fontsize=18)
    
    #change tick label size
    ax.tick_params(axis='x', labelsize=14)
    ax.tick_params(axis='y', labelsize=14)
    
    #add legend and show
    plt.legend(loc='best',fontsize=16)
    title = 'Lattice functions for '+ opts.lattice_name
    plt.title(title, y=1.02, fontsize=20)

    fig1 = plt.gcf()
    plt.show()
    
    if opts.save:
        name = '_'.join(opts.lf_fns) +'_'+opts.lattice_name+ '.pdf'
        if opts.relpath:
            name = os.path.join(opts.relpath,name)
        fig1.savefig(name, bbox_inches='tight')
        

########################################### Getters ###############################################

def get_lf_fns(lattice, lattice_simulator):
    '''
    Return the lattice functions for every element in the lattice, assuming periodicity
    '''
    #define the lattice function names
    lf_names = ("beta_x", "alpha_x", "beta_y", "alpha_y", 
            "psi_x", "psi_y","D_x", "Dprime_x", "D_y", "Dprime_y")
    #construct an empty array of dictionaries corresponding to each lattice element
    lf_info = [{}]
    #loop through lattice elements
    index=0
    for element in lattice.get_elements():
        index += 1
        #get lattice functions for a given element
        lattice_functions = lattice_simulator.get_lattice_functions(element)
        #define dictionary for this element
        lf = {}
        lf['name'] = element.get_name()
        lf['s'] = lattice_functions.arc_length
        lf['length'] = element.get_length()
        #loop through lattice functions for the element
        for lfname in lf_names:
            lf[lfname] = getattr(lattice_functions,lfname)
        #append individual element to array
        lf_info.append(lf)
        
    #construct initial element lf_info[0]
    lf_info[0]['s'] = 0.0
    lf_info[0]['name'] = lattice.get_name()
    lf_info[0]['length'] = 0.0
    lf_info[0]['psi_x'] = 0.0
    lf_info[0]['psi_y'] = 0.0
    #Take lattice functions from the final element to ensure periodicity
    for lfname in lf_names:
        lf_info[0][lfname] = lf_info[-1][lfname]
        
    return lf_info


def get_sliced_lf_fns(lattice, lattice_simulator):
    '''
    Return the lattice functions for every slice in the lattice simulator, assuming periodicity
    '''
    #define the lattice function names
    lf_names = ("beta_x", "alpha_x", "beta_y", "alpha_y", 
            "psi_x", "psi_y","D_x", "Dprime_x", "D_y", "Dprime_y")
    #construct an empty array of dictionaries corresponding to each lattice element
    lf_info = [{}]
    #loop through lattice elements
    index=0
    for aslice in lattice_simulator.get_slices():
        index += 1
        #get lattice functions for a given element
        lattice_functions = lattice_simulator.get_lattice_functions(aslice)
        ele = aslice.get_lattice_element()
        #define dictionary for this element
        lf = {}
        lf['name'] = ele.get_name()
        lf['s'] = lattice_functions.arc_length
        lf['length'] = ele.get_length()
        #loop through lattice functions for the element
        for lfname in lf_names:
            lf[lfname] = getattr(lattice_functions,lfname)
        #append individual element to array
        lf_info.append(lf)
        
    #construct initial element lf_info[0]
    lf_info[0]['s'] = 0.0
    lf_info[0]['name'] = lattice.get_name()
    lf_info[0]['length'] = 0.0
    lf_info[0]['psi_x'] = 0.0
    lf_info[0]['psi_y'] = 0.0
    #Take lattice functions from the final element to ensure periodicity
    for lfname in lf_names:
        lf_info[0][lfname] = lf_info[-1][lfname]
        
    return lf_info


########################################### Main Functions ################################################
        

def plot_sliced_lattice_functions(opts):
    '''Plots the lattice functions for a desired lattice simulator, specified by opts parameters.'''
    
    lf_array = get_sliced_lf_fns(opts.lattice, opts.lattice_simulator)
    
    lf_plot(lf_array, opts)     
    
    

def plot_lattice_functions(opts):
    
    '''Plots the lattice functions for a desired lattice, specified by opts parameters.'''
    
    lf_array = get_lf_fns(opts.lattice, opts.lattice_simulator)
    
    lf_plot(lf_array, opts)   

    
if __name__ == '__main__':
    
    #handle arguments
    lattice = sys.argv[0]
    #define lattice functions being plotted - s is always independent coordinate
    lf_fns = sys.argv[1:]
    
    #get lattice function array
    lf_array = get_lf_fns(lattice)
    lf_plot(lf_array, lf_fns)