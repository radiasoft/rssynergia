import sys, os
import synergia
from mpi4py import MPI
import synergia_workflow
import numpy as np
import rssynergia
from rssynergia.base_diagnostics import workflow


def get_phase_advance(ls, ele0, ele1, x=True):
    '''
    Returns the phase_advance between two elements in a lattice, given a fixed lattice simulator construction.
    Also prints the position of each element and Dispersion in the x-plane at that position.

    Arguments:
        ls (synergia.simulation.lattice_simulator): A Synergia lattice simulator object
        ele0 (str): first element name
        ele1 (str): second element name
        x (Optional[bool]): Whether to look at x phase (True) or y phase (False). Defaults to True.

    Returns:
        advance (float): the phase advance between ele0 and ele1.


    Assumptions:
        - Elements ele0 and ele1 must be unique!
        - Elements ele0 and ele1 must be thin elements.
    '''


    #starting null values
    e0check = False
    e1check = False

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


def make_drift(l):
    """Construct the idealized 6x6 R-matrix for a drift of length l.
    Dispersive components (R-56) are not included.

    Arguments:
        l (float): length of the drift

    Returns:
        drift_R (ndArray): A 6x6 matrix
    """

    drift_R = np.identity(6)
    drift_R[0,1] += l
    drift_R[2,3] += l

    return drift_R

def make_R(f):
    '''Constructs a 6x6 R-matrix focussing in both planes with focal length f.

    Arguments:
        f (float): focal length in each plane.

    Returns:
        m_R (ndArray): A 6x6 matrix
    '''
    mat_R = np.identity(6)
    mat_R[1,0] -= 1/f
    mat_R[3,2] -= 1/f
    mat_R[4,5] += 3.2
    return mat_R

def make_quad(k,l,K=None):
    """Constructs a 6x6 R-matrix for a quadrupole element.

    Arguments:
        k (float): graident of the focusing fields
        l (float): length of the magnet
        K (Optional[float]): integrated strength of the magnet. Overrides k,l inputs. Defaults to None.

    Returns:
        quad_R (ndArray): A 6x6 matrix
    """

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

def get_sextupoles(lattice, nonzero=True):
    '''Return a list of positive and negative sextupoles from a given Synergia lattice.

    Arguments
        lattice (Synergia.latttice.lattice): A Synergia lattice object
        nonzero (Optional[bool]): True if returning only sextupoles with nonzero values. Defaults to True.

    Returns:
    (p_six, n_six) : a length-2 tuple containing lists of positive and negative sectupoles

    '''

    p_six = []
    n_six = []
    last = 'n_six'

    for elem in lattice.get_elements():
        if elem.get_type() == "sextupole":
            k2 = elem.get_double_attribute("k2")
            if k2 > 0.0:
                #p_six.append(elem)
                p_six.append(elem.get_name())
            elif k2 < 0.0:
                #n_six.append(elem)
                n_six.append(elem.get_name())
            elif k2 == 0:
                if nonzero:
                    pass
                else:
                    #sextupole strength is 0, so split elements between positive and negative list
                    if last == 'n_six':
                        #p_six.append(elem)
                        p_six.append(elem.get_name())
                        last = 'p_six'
                    elif last == 'p_six':
                        #n_six.append(elem)
                        n_six.append(elem.get_name())
                        last = 'n_six'

    return (p_six, n_six)

def print_strengths(elemslist, unique=True):
    '''Print the strengths of quads/sextupoles from a list of lattice elements.

    Arguments:
        elemslist(list): a list of Synergia lattice elements (from lattice.get_elements() usually)
        unique (Optional[bool])- If True, print only unique elements. Defaults to True.

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


def get_unique_elements(lattice, autodrifts = False):
    '''
    Return a dictionary of unique elements in the lattice.

    This is essentially a subset of lattice.get_elements(), the dictionary is index-able by element name.

    Arguments:
        lattice (Synergia.lattice.lattice): A synergia lattice object
        autodrifts (Optional[bool]): If True, includes the auto-generated drift elements. Defaults to False.


    Returns:
        subelems_dict (dict): A dictionary of lattice element objects with element names as keys.

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
                if autodrifts == True: #only append auto drifts if flag specified
                    subelems_dict[str(name)] = elem
                    subelems.append(elem)
                    names.append(name)

    return subelems_dict


###################### Chromaticity Correction Scripts ####################################

def chromaticity_adjust(l_s, p_list, n_list, cx, cy = None):
    '''
    Adjust a lattice to meet a desired chromaticity goal in each plane.

    Computes a cost equal to the sum of squared values of sexutpole strengths for the combined
    list of sextupole strengths, prints that cost, and returns it.

    Arguments:
        l_s (Synergia.lattice.lattice): Synergia lattice simulator object
        p_list (list): list of positive sextupole lattice element names which should be varied
        n_list (list): list of negative sextupole lattice element names which should be varied
        cx (float): goal value for horizontal chromaticity
        cy (Optional[float]) - goal value for vertical chromaticity. Default value is cy=cx.

    Returns:
        cost (float): The sum of squared values of sextupole strengths across the list of sextupoles strengths.
    '''

    c_tol = 1.e-3
    c_tol = 1.
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
    #Note - default max steps is 6
    l_s.adjust_chromaticities(cx_goal,cy_goal,p_s,n_s,c_tol,100)

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

########################### Compare lattice functions at entrance##############################
def get_starting_lf(lattice_simulator):
    '''Return the lattice functions at the starting position of the lattice simulator

    Arguments:
        lattice_simulator (Synergia.simulation.lattice_simulator): A Synergia lattice simulator object

    Returns:
        lf.betax: beta functions in the x plane as calculated by the lattice simulator
        lf.betay: beta functions in the y plane as calculated by the lattice simulator
        lf.alphax: alpha functions in the x plane as calculated by the lattice simulator
        lf.alphay: alpha functions in the y plane as calculated by the lattice simulator

    '''


    slices = lattice_simulator.get_slices()
    lf = lattice_simulator.get_lattice_functions(slices[0])

    print "Initial starting lattice functions: betax = {}, betay = {}, alphax = {}, alphay = {}".format(lf.beta_x,lf.beta_y,lf.alpha_x,lf.alpha_y)

    return lf.beta_x,lf.beta_y,lf.alpha_x,lf.alpha_y


###########################Set up stepper and lattice simulator##############################
def generate_stepper(lattice, coll_operator, opts):
    """
    Generates and returns a Synergia stepper given a lattice object and a sequence of options.

    Arguments:
        lattice (Synergia.lattice.lattice): A Synergia lattice object
        coll_operator (Synergia.simulaton.operator): A Synergia collective operator object. May be a dummy operator.
        opts (options.Options): A Synergia options instance

    Returns:
        stepper (Synergia.simulation.stepper): A Synergia stepper object

    """


    requested_stepper = opts.stepper

    if requested_stepper == "splitoperator":
        print "Using split-operator stepper with ", opts.steps, " steps/turn"
        stepper = synergia.simulation.Split_operator_stepper(lattice, opts.map_order, coll_operator, opts.steps)

    elif requested_stepper == "soelements":
        print "Using split-operator stepper elements with ", opts.steps_per_element, " steps/element"
        stepper = synergia.simulation.Split_operator_stepper_elements(lattice, opts.map_order, coll_operator, opts.steps_per_element)

    elif requested_stepper == "independent":
        print "Using independent-operator stepper with ", opts.steps, " steps/turn"
        stepper = synergia.simulation.Independent_stepper(lattice, opts.map_order, opts.steps)

    elif requested_stepper == "elements":
        print "Using step-by-elements-operator stepper with ", opts.steps_per_element, " steps/element"
        stepper = synergia.simulation.Independent_stepper_elements(lattice, opts.map_order, opts.steps_per_element)

    else:
        raise RuntimeError, "stepper %s invalid, must be either 'splitoperator', 'independent' or 'elements'"%requested_stepper


    return stepper


###########################Lattice element string adjustment##############################
def set_lattice_element_type(lattice, opts):
    """
    Adjusts the string attributes of each element in the lattice to reflect the desired extraction type and propagation type.
    May also be used to set an aperture for all elements.


    Arguments:
        lattice (Synergia.lattice.lattice): A Synergia lattice object
        opts (options.Options): A Synergia options instance

    Returns:
        lattice (Synergia.lattice.lattice): The Synergia lattice object, updated with string attributes
    """


    for elem in lattice.get_elements():
        #set an aperture if specified
        if opts.radius:
            elem.set_double_attribute("circular_aperture_radius", opts.radius)

        # set the propagation method.
        # extractor_type may be "chef_maps", "chef_mixed", or chef_propagate
        if opts.use_maps == "taylor":
            pass
        elif opts.use_maps == "all":
            elem.set_string_attribute("extractor_type", "chef_map")
        elif opts.use_maps == "none":
            elem.set_string_attribute("extractor_type", "chef_propagate")
        elif opts.use_maps == "nonrf":
            elem.set_string_attribute("extractor_type", "chef_mixed")
        elif opts.use_maps == "onlyrf":
            if elem.get_type() == "rfcavity":
                elem.set_string_attribute("extractor_type", "chef_map")
            else:
                elem.set_string_attribute("extractor_type", "chef_propagate")
        else:
            raise RuntimeError, "bad options for use_maps: %d"%opts.use_maps

    return lattice


#############################Lattice Comparison Script ###################################
def compare_lattices(lattice1, lattice2):

    '''
    Compares the elements between two lattices.

    Lists elements which are unique to each lattice, and further lists discrepancies in elements
    that share the same name. Returns no information for lattice elements which are the same.

    Arguments:
        lattice1 (Synergia.lattice.lattice): A Synergia lattice object
        lattice2 (Synergia.lattice.lattice): A Synergia lattice object

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
    '''
    Return the 6x6 R-matrices for each step in a default Synergia stepper for a given lattice.

    This function runs through every step in the stepper, extracts the operators for each step,
    and processes these using Synergia's fast mapping operation to obtain 6x6 matrices for each operator.
    A list comprising all such maps is then returned for the input lattice object.

    Arguments:
        lattice (Synergia.lattice.lattice): A Synergia lattice object

    Returns:
        maps (list): A list of 6x6 matrices obtained from Synergia's fast mapping operation for each step.

    '''

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
    """
    Multiplies together a list of maps to construct a combined map.

    Arguments:
        maplist (list): A list of 6x6 arrays representing R-matrices.

    Returns:
        final (ndArray): The product of the matrices in maplist.
    """
    copylist = list(maplist) #copy the list to preserve the original
    dim = copylist[0].shape[0]
    final = np.identity(dim)
    while len(copylist) > 0:
        mult = copylist.pop()
        final = np.dot(mult,final)
    return final
