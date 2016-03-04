import sys, os
import synergia
from mpi4py import MPI
import synergia_workflow
import numpy as np


def get_phase_advance(ls, ele0, ele1, x=True):
    '''
    Returns the phase_advance between two elements in a lattice, given a fixed lattice simulator construction.
    
    Arguments:
        -ls - lattice simulator
        -ele0 - first element name 
        -ele1 - second element name
        
    Assumptions:
        - Elements ele0 and ele1 must be unique!
        - Elements ele0 and ele1 must be thin elements.
    '''
    
    
    #starting null values
    e0check = False
    e1check = False
    
    #for ele in ls.get_lattice().get_elements():
    #    if ele.get_name() == ele0:
    #        e0 = ele
    #    if ele.get_name() == ele1:
    #        e1 = ele
            
    #if not e0:
    #    print "Cannot find element: {}".format(ele0)
    #    return 0
    #elif not e1:
    #    print "Cannot find element: {}".format(ele1)
    #    return 0
            
    parts = ls.get_slices()
    arclengths = [ls.get_lattice_functions(part).arc_length for part in parts]
    
    for index,part in enumerate(parts):
        if part.get_lattice_element().get_name() == ele0:
            #print "Found {}".format(ele0)
            e0check = True
            lf = ls.get_lattice_functions(part)
            if x == True:
                begin = lf.psi_x
            else:
                begin = lf.psi_y
                
            endpoint = arclengths[index] #endpoint of slice
            print "Element {} is at position {} with dispersion Dx = {}".format(ele0,endpoint,lf.D_x)
            
        if part.get_lattice_element().get_name() == ele1:
            #print "Found {}".format(ele1)
            e1check = True
            lf = ls.get_lattice_functions(part)
            if x == True:
                end = lf.psi_x
            else:
                end = lf.psi_y
                
            endpoint = arclengths[index] #endpoint of slice
            print "Element {} is at position {} with dispersion Dx = {}".format(ele1,endpoint,lf.D_x)
            
    if not e0check:
        print "Element {} not found".format(ele0)
        return 0
    elif not e1check:
        print "Element {} not found".format(ele1)
        return 0
            
    advance = (end - begin)/np.piIOTA
    #advance = 1
    
    if x:
        coord = 'x'
    else:
        coord = 'y'
    
    if advance < 0:
        print "Elements specified out of order"
        print "Elements {} and {} are separated by a phase advance in {} of {}".format(ele0, ele1,coord, advance)
        return -1.*advance
        
    print "Elements {} and {} are separated by an phase advance in {} of {}".format(ele0, ele1,coord, advance)
    
    return advance
    

def make_drift(length):

    drift_R = np.identity(6)
    drift_R[0,1] += length
    drift_R[2,3] += length
    #Include the 5-6 component for completeness
    
    return drift_R

def make_R(f):
    '''Make an R-matrix focussing in both planes to simulate the IOTA lattice outside of the NLL elements'''
    mat_R = np.identity(6)
    mat_R[1,0] -= 1/f
    mat_R[3,2] -= 1/f
    mat_R[4,5] += 3.2
    return mat_R

def make_quad(k,l,K=None):
    
    if K:
        quad_R = np.identity(6)
        quad_R[1,0] -= K
        quad_R[3,2] += K

    else:
        quad_R = np.identity(6)
        quad_R[1,0] -= k*l
        quad_R[3,2] += k*l
        
    
    return quad_R  

def make_chef(lattice, lattice_simulator=None):
    '''Make all elements of a lattice chef-propogated and update the lattice simulator.'''
    
    elems = lattice.get_elements()
    #nlelems = []
    for elem in elems:
        if elem.get_type() == 'nllens':
            #nlelems.append(elem)
            elem.set_string_attribute("extractor_type", "chef_propagate")
        else:
            elem.set_string_attribute("extractor_type", "chef_propagate")
    if lattice_simulator:
        lattice_simulator.update()


def get_fd_quads(lattice):
    '''Return a list of focussing and defocussing quads in a lattice
    
    Arguments
        -lattice - Synergia lattice object
        
    Outputs a pair of lists -  ( [list of focussing quad elements], [list of defocussing quad elements] )
    
    '''
    
    f_quads = []
    d_quads = []
    for elem in lattice.get_elements():
        if elem.get_type() == "quadrupole":
            k1 = elem.get_double_attribute("k1")
            if k1 > 0.0:
                f_quads.append(elem)
            elif k1 < 0.0:
                d_quads.append(elem)
    return (f_quads, d_quads)
    
def get_sextupoles(lattice):
    '''Return a list of positive and negative sextupoles
    
    Arguments
        -lattice - Synergia lattice object
        
    Outputs a pair of lists -  ( [list of positive sextupoles], [list of negative sextupoles] )
    
    '''
    
    p_six = []
    n_six = []
    last = 'n_six'
    
    for elem in lattice.get_elements():
        if elem.get_type() == "sextupole":
            k2 = elem.get_double_attribute("k2")
            if k2 > 0.0:
                p_six.append(elem)
            elif k2 < 0.0:
                n_six.append(elem)
            elif k2 == 0:
                'sextupole strength is 0, so split elements between positive and negative list'
                if last == 'n_six':
                    p_six.append(elem)
                    last = 'p_six'
                elif last == 'p_six':
                    n_six.append(elem)
                    last = 'n_six'
 
    return (p_six, n_six)


#quick helper method
def print_strengths(elemslist, unique=True):
    '''Print the strengths of quads/sextupoles from a list of lattice elements
    
    Arguments:
        elemslist - the list of elements (from lattice.get_elements() usually)
        
    Optional:
        unique - (defaults to True) print only unique elements
    
    '''
    strengths = []
    for elem in elemslist:
        if elem.get_type() == "quadrupole":
            strength = elem.get_double_attribute("k1")
            if unique:
                #only print new strengths
                if strength not in strengths:
                    print elem.get_name() + ' K: ' + str(elem.get_double_attribute("k1"))
                    strengths.append(strength)
                else:
                    pass
            else:
                #print all strengths
                print elem.get_name() + ' K: ' + str(elem.get_double_attribute("k1"))
        #do the same for sextupoles
        if elem.get_type() == "sextupole":
            strength = elem.get_double_attribute("k2")
            if unique:
                #only print new strengths
                if strength not in strengths:
                    print elem.get_name() + ' K2: ' + str(elem.get_double_attribute("k2"))
                    strengths.append(strength)
                else:
                    pass
            else:
                #print all strengths
                print elem.get_name() + ' K2: ' + str(elem.get_double_attribute("k2"))             


def get_unique_elements(lattice, autodrifts = None):
    '''
    Return a dictionary of unique elements in the lattice
    
    This will be a subset of lattice.get_elements(), the dictionary is index-able by element name.
    
    Arguments:
    
    -lattice - the synergia lattice object
    
    -autodrifts (default None) - specify this flag to include the drifts created by 
    synergia to automatically bridge elements in a sequence.
    
    '''
    
    elems = lattice.get_elements()
    
    subelems = []
    subelems_dict = {}
    enames = []
    for elem in elems:
        name = elem.get_name()
        if name in enames:
            pass
        else:
            if not name.startswith('auto_drift'):
                subelems_dict[str(name)] = elem
                subelems.append(elem)
                enames.append(name)
            else:
                if autodrifts: #only append auto drifts if flag specified
                    subelems_dict[str(name)] = elem
                    subelems.append(elem)
                    names.append(name)    
            
    return subelems_dict
    

###################### Chromaticity Correction Script ####################################
from base_diagnostics import workflow

def DEPRECATED_adjust_chromaticity(lattice, maglist,):

    '''This guy is very much a work in progress!'''

    #Get the uncorrected lattice!
    adjusted_lattice_simulator = stepper2.get_lattice_simulator()
    adjusted_lattice = adjusted_lattice_simulator.get_lattice()

    cy_orig = adjusted_lattice_simulator.get_vertical_chromaticity()
    cx_orig = adjusted_lattice_simulator.get_horizontal_chromaticity()

    cx_goal = cy_orig
    cy_goal = cy_orig

    c_tol = 1.0e-5; #need to define a tolerance for the fit


    #grab the sextupoles which are being used (nonzero) in the chromaticity corrected lattice
    p_six_c, n_six_c = latticework.get_sextupoles(lattice)
    p_six_use = [ele.get_name() for ele in p_six_c if ele.get_double_attribute("k2") > 0]
    n_six_use = [ele.get_name() for ele in n_six_c if ele.get_double_attribute("k2") < 0]

    #print the original strengths of the things that will be adjusted
    #print [": ".join([ele.get_name(), str(ele.get_double_attribute("k2"))]) for ele in p_six_c]
    #print [": ".join([ele.get_name(), str(ele.get_double_attribute("k2"))]) for ele in n_six_c]

    #grab the sextupoles from the adjusted lattice
    p_six, n_six = latticework.get_sextupoles(adjusted_lattice)

    #filter to make sure we only keep the ones that we want to adjust
    p_s = [ele for ele in adjusted_lattice.get_elements() if ele.get_name() in p_six_use]
    n_s = [ele for ele in adjusted_lattice.get_elements() if ele.get_name() in n_six_use]

    #use these magnets to adjust the chromaticity
    adjusted_lattice_simulator.adjust_chromaticities(cx_goal,cy_goal,p_s,n_s,c_tol)

    #update the lattice
    adjusted_lattice_simulator.update()
    adjusted_lattice = adjusted_lattice_simulator.get_lattice()

    #Force define original lattice simulator
    lattice2 = synergia.lattice.MadX_reader().get_lattice("iota", "/Users/ncook/Synergia_Tests/lattices/Iota6-6/lattice_2IO.madx")
    stepper2 = synergia.simulation.Independent_stepper_elements(lattice2, opts.map_order, opts.steps_per_element)
    lattice_simulator2 = stepper2.get_lattice_simulator()


    #Print stuff for original lattice
    print "Original horizontal chromaticity: {}".format(lattice_simulator2.get_horizontal_chromaticity())
    print "Original vertical chromaticity: {}".format(lattice_simulator2.get_vertical_chromaticity())

    #print the original strengths of the things that will be adjusted
    p_six_orig, n_six_orig = latticework.get_sextupoles(lattice2)
    latticework.print_strengths(p_six_orig, False)
    latticework.print_strengths(n_six_orig, False)

    print ""
    #Print stuff for new re-tuned lattice
    print "New horizontal chromaticity: {}".format(adjusted_lattice_simulator.get_horizontal_chromaticity())
    print "New vertical chromaticity: {}".format(adjusted_lattice_simulator.get_vertical_chromaticity())
    p_six_new, n_six_new = latticework.get_sextupoles(adjusted_lattice)
    latticework.print_strengths(p_six_new, False)
    latticework.print_strengths(n_six_new, False)
    #print [": ".join([ele.get_name(), str(ele.get_double_attribute("k2"))]) for ele in p_six_new]
    #print [": ".join([ele.get_name(), str(ele.get_double_attribute("k2"))]) for ele in n_six_new]




def chromaticity_adjust(l_s, p_list, n_list, cx, cy = None):
    '''
        -lattice - synergia lattice
        -p_list - list of positive sextupole lattice element names which should be varied
        -n_list - list of negative sextupole lattice element names which should be varied 
        -C0 - a value for C_x, C_y we want to reach
        
        DEPRECTATED:
            -p_s - list of positive-valued sextupoles to adjust
            -n_s - list of negative-values sextupoles to adjust
    '''
    
    c_tol = 1.e-3
    cx_goal = cx
    if cy:
        cy_goal = cy
    else:
        cy_goal = cx_goal
    
    #make sure we grab the correct elements to adjust
    lattice = l_s.get_lattice()
    p_s = [ele for ele in lattice.get_elements() if ele.get_name() in p_list]
    n_s = [ele for ele in lattice.get_elements() if ele.get_name() in n_list]
    
    #adjust the chromaticity
    l_s.adjust_chromaticities(cx_goal,cy_goal,p_s,n_s,c_tol)
    
    cv_final = l_s.get_vertical_chromaticity()
    cx_final = l_s.get_horizontal_chromaticity()
    print 'Adjusted chromaticities to new values: Cx = {} , Cy = {}'.format(cx_final, cv_final)
    print ''
    
    for ele in p_s:
        print "New k2 value for element {}: {}".format(ele.get_name(), ele.get_double_attribute('k2'))
    for ele in n_s:
        print "New k2 value for element {}: {}".format(ele.get_name(), ele.get_double_attribute('k2'))
    
    #calculate Cost function
    cost = 0
    cost += np.sum([np.abs(ele.get_double_attribute('k2')) for ele in p_s]) #add p_s contributions
    cost += np.sum([np.abs(ele.get_double_attribute('k2')) for ele in n_s]) #add n_s contributions
    
    print 'Calculated Cost Function: {}'.format(cost)
    
    return cost




#############################Lattice Comparison Script ###################################
def compare_lattices(lattice1, lattice2):

    '''
    
    Compares the elements in a lattice. Lists elements which are unique to each lattice, and further lists discrepancies
    in elements that share the same name. Returns no information for lattice elements which are the same.
    
    
    '''


    elems1 = get_unique_elements(lattice1)
    elems2 = get_unique_elements(lattice2)
    l1 = len(elems1)
    l2 = len(elems2)
    
    print "Lattice 1 has " + str(l1) + " elements."
    print "Lattice 2 has " + str(l2) + " elements."
    print
    
    l1_keys = elems1.keys()
    l2_keys = elems2.keys()
    
    l_shared = set(l1_keys).intersection(l2_keys) #return shared items for each lattice
    l1_unique = set(l1_keys).difference(l2_keys) #return items unique to lattice 1
    l2_unique = set(l2_keys).difference(l1_keys) #return items unique to lattice 2
    
    print "Lattice 1 has " + str(len(l1_unique)) + " unique elements :"
    print [element for element in l1_unique]
    print
    print "Lattice 2 has " + str(len(l2_unique)) + " unique elements :"
    print [element for element in l2_unique]
    print    
    print "The lattices have " + str(len(l_shared)) + " shared elements."
    print
    
    counter = 0
    #for index, elem in enumerate(elems1):
    #for key, value in elems1.iteritems():
    #now loop through all of the share elements
    for key in l_shared:
        elem = elems1[key]
        elem2 = elems2[key]
        if elem.get_type() == elem2.get_type(): #check type
            da1 = elem.get_double_attributes()
            da2 = elem2.get_double_attributes()
            
            da1keys = da1.keys()
            da2keys = da2.keys()
            
            if da1keys == da2keys: #check double attribute keys
                for dkey in da1keys: #check double attribute values
                    if not da1[dkey] == da2[dkey]:
                        print 'Variation in element double attribute for element ' + str(dkey) + ' of type ' + elem.get_type()
                        print 'Lattice 1 element ' + str(elem.get_name()) + ' has value ' + str(da1[dkey])
                        print 'Lattice 2 element ' + str(elem2.get_name()) + '  has value ' + str(da2[dkey])
                        counter += 1
            else:
                #ignore the 'deposited_charge' attribute
                #print set(da1keys).symmetric_difference(da2keys)
                if not set(da1keys).symmetric_difference(da2keys) == set(['deposited_charge']):
                    print 'Variation in element double attribute keys for element ' + str(key) + ' of type ' + elem.get_type()
                    print 'Lattice 1 element ' + str(elem.get_name()) + ' has double attributes ' + str([key for key in da1.keys()])
                    print 'Lattice 2 element ' + str(elem2.get_name()) + '  has double attribues ' + str([key for key in da2.keys()])
                    counter += 1
                else:
                    pass 
        else:
            print 'Variation in element type for element ' + str(key)
            print 'Lattice 1 element ' + str(elem.get_name()) + ' has type ' + elem.get_type()
            print 'Lattice 2 element ' + str(elem2.get_name()) + ' has type ' + elem2.get_type()
            counter += 1
            
    
    print 'There are ' + str(counter) + ' element specific variations in the shared elements of the lattices.'




def get_maps(lattice):
    '''Return the matrices for each step in the stepper'''
    
    for element in lattice.get_elements():
        element.set_string_attribute("extractor_type", "chef_map")
    map_order = 1
    lattice_simulator = synergia.simulation.Lattice_simulator(lattice,
                                                              map_order)
    steps_per_element = 1
    stepper = synergia.simulation.Independent_stepper_elements(lattice_simulator,
                                                                steps_per_element)
    
    maps = []
    
    stepper.force_update_operations_no_collective()
    for step in stepper.get_steps():
        for operator in step.get_operators():
            if operator.get_type() == 'independent':
                io = synergia.simulation.as_independent_operator(operator)
                slices = io.get_slices()
                #print "%s:" % slices[0].get_lattice_element().get_name()
                if len(slices) > 1:
                    raise RuntimeError, "found multiple slices in an independent operator"
                for operation in io.get_operations():
                    if operation.get_type() == 'fast_mapping':
                        fmo = synergia.simulation.as_fast_mapping_operation(operation)
                        dense_mapping = synergia.simulation.Dense_mapping(fmo.get_fast_mapping())
                        #print "order 0:"
                        #print dense_mapping.get_constant_term()
                        #print "order 1:"
                        #print dense_mapping.get_linear_term()
                        #etype = io.get_slices()[0].get_lattice_element().get_type()
                        #ename = io.get_slices()[0].get_lattice_element().get_name()
                        dmap = dense_mapping.get_linear_term()
                        dmap_copy = np.copy(dmap)
                        #dmmap = dense_mapping.get_linear_term_mad()                        
                        maps.append(np.asarray(dmap_copy))
                        #if etype == 'quadrupole':
                            #dmap = dense_mapping.get_linear_term()
                            #dmap_copy = np.copy(dmap)
                            #dmmap = dense_mapping.get_linear_term_mad()
                            #print ename + ":"
                            #print dmap
                            #if dmap.all() == dmmap.all():
                                #print "mad map is same"
                            #maps.append(np.asarray(dmap_copy))
    return maps
    
def combine_maps(maplist):
    copylist = list(maplist) #copy the list to preserve the original
    dim = copylist[0].shape[0]
    final = np.identity(dim)
    while len(copylist) > 0:
        mult = copylist.pop()
        final = np.dot(mult,final)
    return final
